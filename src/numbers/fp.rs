use std::ops::*;

use crate::numbers::*;
use libbgs_util::*;

/// An integer modulo `P`.
///
/// # Example
/// ```
/// use libbgs::numbers::FpNum;
/// let x = FpNum::<7>::from(5);
/// let y = FpNum::<7>::from(6);
/// let z = x * y;
/// assert_eq!(z, FpNum::from(30 % 7));
/// ```
///
/// For more informtion on the internal representation, see: Montgomery, Peter (April 1985). "Modular Multiplication Without Trial
/// Division". Mathematics of Computation. 44 (170): 519-521.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FpNum<const P: u128>(u128);

impl<const P: u128> FpNum<P> {
    /// The constant 0.
    pub const ZERO: FpNum<P> = FpNum::from_u128(0);

    // Once const trait impls are stabalized, this can be replaced with a call to the pow method
    // from the GroupElem trait.
    // Until then, we copy + paste the code from there... not very DRY of me.
    /// The constant $2^{-1}$.
    pub const TWO_INV: FpNum<P> = FpNum::from_u128(2).const_pow(P - 2);

    /// Returns the Legendre symbol of `a` modulo `P`, i.e.,
    /// $$\left(\frac{a}{p}\right)_L = a^{\frac{p - 1}{2}} \mod p$$.
    pub const fn legendre(&self) -> FpNum<P> {
        self.const_pow((P - 1) / 2)
    }

    /// Calculates this number's square root, if it is a quadratic residue; otherwise, returns
    /// `None`.
    pub const fn int_sqrt(&self) -> Option<FpNum<P>> {
        if self.0 == FpNum::<P>::from_u128(0).0 {
            return Some(FpNum::ZERO);
        }

        // Tonelli algorithm
        let mut s = 0;
        let mut q = P - 1;
        while q % 2 == 0 {
            s += 1;
            q /= 2;
        }

        let mut i = 1;
        let z = loop {
            let z = FpNum::from_u128(standard_affine_shift(P, i));
            if z.legendre().0 == FpNum::<P>::from_u128(P - 1).0 {
                break z;
            }
            i += 1;
        };
        let mut c = z.const_pow(q);
        let mut r = self.const_pow((q + 1) / 2);
        let mut t = self.const_pow(q);
        let mut m = s;

        while t.0 != FpNum::<P>::ONE.0 {
            let mut temp = t;
            let mut i = 0;
            while temp.0 != FpNum::<P>::ONE.0 {
                temp = temp.const_mul(&temp);
                i += 1;
            }
            if i == m {
                return None;
            }
            let b = c.const_pow(1 << (m - i - 1));

            r = r.const_mul(&b);
            c = b.const_mul(&b);
            t = t.const_mul(&c);
            m = i;
        }
        Some(r)
    }

    /// Returns a quadratic nonresidue modulo `p`.
    pub const fn find_nonresidue() -> FpNum<P> {
        if P % 4 == 3 {
            FpNum::<P>::from_u128(P - 1)
        } else if P % 8 == 3 || P % 8 == 5 {
            FpNum::<P>::from_u128(2)
        } else {
            let mut i = 0;
            let res = loop {
                let a = standard_affine_shift(P, i);
                if intpow::<P>(a, (P - 1) / 2) == P - 1 {
                    break a;
                }
                i += 1;
            };
            FpNum::<P>::from_u128(res)
        }
    }

    /// Returns the Montgomery representation of this number.
    #[inline(always)]
    pub const fn raw(&self) -> u128 {
        self.0
    }


    /// Converts a `u128` into its Montgomery representation.
    /// This operation is expensive.
    pub(crate) const fn from_u128(src: u128) -> FpNum<P> {
        let r2 = long_multiply::<P>(Self::R, Self::R);
        FpNum::<P>::redc2(carrying_mul(src, r2))
    }

    // Elements are represented as aR + N in this system.
    // Chosen to be a power of two to turn multiplication / division into shift operations.
    const R: u128 = {
        let mut n = P.next_power_of_two();
        while n.ilog2() % 4 != 0 {
            n <<= 1;
        }
        n
    };

    // Mask used for taking elements modulo R. Since R is a power of 2,
    // x % R == x & MASK
    const MASK: u128 = Self::R - 1;

    // The length of R
    const SHIFT: u32 = Self::R.ilog2();

    // "Magic" number R' such that RR' is congruent to 1 mod N
    const MAGIC: u128 = {
        let mut r0 = Self::R;
        let mut r1 = P;
        let mut t0: i128 = 0;
        let mut t1: i128 = 1;
        let mut s0: i128 = 1;
        let mut s1: i128 = 0;
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (t0, t1) = (t1, t0 - (q as i128) * t1);
            (s0, s1) = (s1, s0 - (q as i128) * s1);
        }
        if t0 < 0 {
            (-t0) as u128
        } else {
            Self::R - (t0 as u128)
        }
    };

    #[allow(dead_code)]
    const fn redc(x: u128) -> FpNum<P> {
        let m = ((x & Self::MASK) * Self::MAGIC) & Self::MASK;
        let t = (x + m * P) >> Self::SHIFT;
        if t >= P {
            FpNum(t - P)
        } else {
            FpNum(t)
        }
    }

    const fn redc2((hi, lo): (u128, u128)) -> FpNum<P> {
        let (_, m) = carrying_mul(lo & Self::MASK, Self::MAGIC);
        let m = m & Self::MASK;
        let (c1, t) = carrying_mul(m, P);
        let (c2, t) = carrying_add(t, lo);
        let t = shrd(t, c1 + c2 + hi, Self::SHIFT as usize);
        if t >= P {
            FpNum(t - P)
        } else {
            FpNum(t)
        }
    }

    // These next several methods are constant implementations of the usual std::ops methods.
    // Once const impls lands in Rust, they will be removed, and the impl std::ops will be made
    // const.
    pub(crate) const fn const_mul(&self, rhs: &FpNum<P>) -> FpNum<P> {
        let (hi, lo) = carrying_mul(self.0, rhs.0);
        FpNum::<P>::redc2((hi, lo))
    }

    pub(crate) const fn const_pow(self, mut n: u128) -> FpNum<P> {
        let mut x = self;
        let mut y = FpNum::<P>::from_u128(1);
        if n == 0 {
            return y;
        }
        while n > 1 {
            if n % 2 == 1 {
                y = y.const_mul(&x);
            }
            x = x.const_mul(&x);
            n >>= 1;
        }
        x.const_mul(&y)
    }
}

impl<S, const P: u128> SylowDecomposable<S> for FpNum<P>
where
    FpNum<P>: Factor<S>,
{
    fn find_sylow_generator(i: usize) -> FpNum<P> {
        match Self::FACTORS[i] {
            (2, 1) => FpNum::from(FpNum::<P>::SIZE),
            (p, t) => (1..FpNum::<P>::SIZE)
                .map(|j| FpNum::from(standard_affine_shift(P, j)))
                .find_map(|c| <FpNum<P> as SylowDecomposable<S>>::is_sylow_generator(&c, (p, t)))
                .unwrap(),
        }
    }
}

impl<const P: u128> GroupElem for FpNum<P> {
    const ONE: Self = FpNum::from_u128(1);
    const SIZE: u128 = P - 1;

    fn multiply(&self, other: &FpNum<P>) -> FpNum<P> {
        *self * *other
    }
}

impl<const P: u128> From<u128> for FpNum<P> {
    fn from(src: u128) -> FpNum<P> {
        FpNum::<P>::from_u128(src)
    }
}

impl<const P: u128> From<FpNum<P>> for u128 {
    fn from(src: FpNum<P>) -> u128 {
        FpNum::<P>::redc2((0, src.0)).0
    }
}

impl<const P: u128> Add<FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, rhs: FpNum<P>) -> FpNum<P> {
        let sum = self.0 + rhs.0;
        if sum >= P {
            FpNum(sum - P)
        } else {
            FpNum(sum)
        }
    }
}

impl<const P: u128> Add<FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, rhs: FpNum<P>) -> FpNum<P> {
        let sum = self.0 + rhs.0;
        if sum >= P {
            FpNum(sum - P)
        } else {
            FpNum(sum)
        }
    }
}

impl<const P: u128> Add<&FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, rhs: &FpNum<P>) -> FpNum<P> {
        let sum = self.0 + rhs.0;
        if sum >= P {
            FpNum(sum - P)
        } else {
            FpNum(sum)
        }
    }
}

impl<const P: u128> Add<&FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, rhs: &FpNum<P>) -> FpNum<P> {
        let sum = self.0 + rhs.0;
        if sum >= P {
            FpNum(sum - P)
        } else {
            FpNum(sum)
        }
    }
}

impl<const P: u128> Sub<FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, rhs: FpNum<P>) -> FpNum<P> {
        if self.0 >= rhs.0 {
            FpNum(self.0 - rhs.0)
        } else {
            FpNum(P + self.0 - rhs.0)
        }
    }
}

impl<const P: u128> Sub<FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, rhs: FpNum<P>) -> FpNum<P> {
        if self.0 >= rhs.0 {
            FpNum(self.0 - rhs.0)
        } else {
            FpNum(P + self.0 - rhs.0)
        }
    }
}

impl<const P: u128> Sub<&FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, rhs: &FpNum<P>) -> FpNum<P> {
        if self.0 >= rhs.0 {
            FpNum(self.0 - rhs.0)
        } else {
            FpNum(P + self.0 - rhs.0)
        }
    }
}

impl<const P: u128> Sub<&FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, rhs: &FpNum<P>) -> FpNum<P> {
        if self.0 >= rhs.0 {
            FpNum(self.0 - rhs.0)
        } else {
            FpNum(P + self.0 - rhs.0)
        }
    }
}

impl<const P: u128> Mul<FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, rhs: FpNum<P>) -> FpNum<P> {
        self.const_mul(&rhs)
    }
}

impl<const P: u128> Mul<FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, rhs: FpNum<P>) -> FpNum<P> {
        self.const_mul(&rhs)
    }
}

impl<const P: u128> Mul<&FpNum<P>> for FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, rhs: &FpNum<P>) -> FpNum<P> {
        self.const_mul(rhs)
    }
}

impl<const P: u128> Mul<&FpNum<P>> for &FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, rhs: &FpNum<P>) -> FpNum<P> {
        self.const_mul(rhs)
    }
}

impl<const P: u128> Mul<FpNum<P>> for u128 {
    type Output = FpNum<P>;
    fn mul(self, rhs: FpNum<P>) -> FpNum<P> {
        rhs.mul(FpNum::from(self))
    }
}

impl<const P: u128> MulAssign<FpNum<P>> for FpNum<P> {
    fn mul_assign(&mut self, rhs: FpNum<P>) {
        *self = self.const_mul(&rhs);
    }
}

impl<const P: u128> Neg for FpNum<P> {
    type Output = FpNum<P>;
    fn neg(self) -> FpNum<P> {
        if self.0 == 0 {
            FpNum(0)
        } else {
            FpNum(P - self.0)
        }
    }
}

impl<const P: u128> Neg for &FpNum<P> {
    type Output = FpNum<P>;
    fn neg(self) -> FpNum<P> {
        if self.0 == 0 {
            FpNum(0)
        } else {
            FpNum(P - self.0)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    impl Factor<Phantom> for FpNum<13> {
        const FACTORS: Factorization = Factorization::new(&[(2, 2), (3, 1)]);
    }

    impl Factor<Phantom> for FpNum<29> {
        const FACTORS: Factorization = Factorization::new(&[(2, 2), (7, 1)]);
    }

    impl Factor<Phantom> for FpNum<BIG_P> {
        const FACTORS: Factorization = Factorization::new(&[
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
    }

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    #[test]
    fn multiplies() {
        let mut x = FpNum::<7>::from(3);
        let five = FpNum::<7>::from(5);
        x = x.multiply(&five);
        assert_eq!(FpNum::ONE, x);
    }

    #[test]
    fn squares() {
        let mut x = FpNum::<7>::from(3);
        x *= x;
        assert_eq!(FpNum::from(2), x);

        let mut x = FpNum::<7>::from(2);
        x *= x;
        assert_eq!(FpNum::from(4), x);
    }

    #[test]
    fn powers_up() {
        let mut x = FpNum::<7>::from(2);
        x = x.pow(5);
        assert_eq!(FpNum::from(4), x);

        let mut x = FpNum::<7>::from(3);
        x = x.pow(3);
        assert_eq!(FpNum::from(6), x);

        let mut x = FpNum::<7>::from(5);
        x = x.pow(6);
        assert!(x == FpNum::ONE);
    }

    #[test]
    fn powers_up_big() {
        let mut x = FpNum::<BIG_P>::from(3);
        x = x.pow(BIG_P - 1);
        assert!(x == FpNum::ONE);
    }

    #[test]
    fn sylow_one_is_one() {
        let one = SylowElem::<Phantom, 2, FpNum<13>>::ONE;
        assert!(one == SylowElem::ONE);
    }

    #[test]
    fn sylow_finds_generators() {
        let g = SylowDecomp::<Phantom, 2, FpNum<29>>::new();
        for i in 0..2 {
            let gen = g.generator(i);
            let d = SylowElem::<Phantom, 2, FpNum<29>>::FACTORS.factor(i);
            test_is_generator_small::<Phantom, 2, FpNum<29>>(gen, d as usize);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let g = SylowDecomp::<Phantom, 7, FpNum<BIG_P>>::new();
        for i in 0..7 {
            let gen = g.generator(i);
            let d = SylowElem::<Phantom, 7, FpNum<BIG_P>>::FACTORS[i];
            test_is_generator_big::<Phantom, 7, FpNum<BIG_P>>(gen, d);
        }
    }

    #[test]
    fn sylow_order() {
        for i in 1..13 {
            let mut x = SylowElem::<Phantom, 2, FpNum<13>>::new([i % 4, i % 3]);
            x = x.pow(x.order());
            assert!(x == SylowElem::ONE);
        }
    }

    #[test]
    fn sylow_order_big() {
        let n = 123456789;
        let mut x = SylowElem::<Phantom, 7, FpNum<BIG_P>>::new(
            SylowElem::<Phantom, 7, FpNum<BIG_P>>::FACTORS
                .factors()
                .iter()
                .map(|(p, d)| n % intpow::<0>(*p, *d as u128))
                .collect::<Vec<u128>>()
                .try_into()
                .unwrap(),
        );
        let or = x.order();
        x = x.pow(or);
        assert!(x == SylowElem::ONE);
    }

    #[test]
    fn calculates_square_roots() {
        let mut nonresidues = 0;
        for x in (1..13).map(|i| FpNum::<13>::from(i)) {
            match x.int_sqrt() {
                None => {
                    nonresidues += 1;
                }
                Some(mut y) => {
                    y = y.pow(2);
                    assert_eq!(x, y);
                }
            }
        }
        assert_eq!(nonresidues, 6);
    }

    #[test]
    fn inverses() {
        for i in 2..13 {
            let mut x = FpNum::<13>::from(i);
            let y = x.clone();
            x = x.inverse();
            assert!(x != FpNum::ONE);
            x = x.multiply(&y);
            assert!(x == FpNum::ONE);
        }
    }

    #[test]
    fn finds_magic() {
        assert_eq!(FpNum::<13>::MAGIC, 11);
    }

    #[test]
    fn test_multiply() {
        let five = FpNum::<13>::from(5);
        let six = FpNum::<13>::from(6);
        let thirty = &five * &six;
        assert_eq!(u128::from(thirty), 30 % 13);
    }

    #[test]
    fn test_multiply_big() {
        let a = 0xABCDEF0123456789;
        let b = 0x9876543210ABCDEF;
        let am = FpNum::<BIG_P>::from(a);
        let bm = FpNum::<BIG_P>::from(b);
        let res = am * bm;
        assert_eq!(u128::from(res), (a * b) % BIG_P);
    }

    #[test]
    fn test_pow() {
        let a = 0xABCDEF;
        let n = 0xFEDCBA;
        let res = FpNum::<BIG_P>::from(a).const_pow(n);
        assert_eq!(u128::from(res), intpow::<BIG_P>(a, n));
    }
}
