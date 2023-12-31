use std::ops::*;

use crate::numbers::*;
use libbgs_util::*;

/// An integer modulo `P`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FpNum<const P: u128>(pub Montgomery<P>);

impl<const P: u128> FpNum<P> {
    /// The constant 0.
    pub const ZERO: FpNum<P> = FpNum(Montgomery::from_u128(0));
    
    // Once const trait impls are stabalized, this can be replaced with a call to the pow method
    // from the GroupElem trait.
    // Until then, we copy + paste the code from there... not very DRY of me.
    /// The constant $2^{-1}$.
    pub const TWO_INV: FpNum<P> = FpNum(Montgomery::from_u128(2).const_pow(P - 2));

    /// Returns the Legendre symbol of `a` modulo `P`, i.e.,
    /// $$\left(\frac{a}{p}\right)_L = a^{\frac{p - 1}{2}} \mod p$$.
    pub fn legendre(&self) -> FpNum<P> {
        self.pow((P - 1) / 2)
    }

    /// Calculates this number's square root, if it is a quadratic residue; otherwise, returns
    /// `None`.
    pub fn int_sqrt(&self) -> Option<FpNum<P>> {
        if self.0 == Montgomery::<P>::from_u128(0) {
            return Some(FpNum::from(0));
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
            let z = FpNum::from(standard_affine_shift(P, i));
            if z.legendre() == (P - 1) {
                break z;
            }
            i += 1;
        };
        let mut c = z.pow(q);
        let mut r = self.pow((q + 1) / 2);
        let mut t = self.pow(q);
        let mut m = s;

        loop {
            if t == 1 {
                return Some(r);
            }
            let mut temp = t;
            let mut i = 0;
            while temp != 1 {
                temp *= temp;
                i += 1;
            }
            if i == m {
                return None;
            }
            let b = c.pow(1 << (m - i - 1));

            r *= b;
            c = b * b;
            t *= c;
            m = i;
        }
    }

    /// Returns a quadratic nonresidue modulo `p`.
    pub const fn find_nonresidue() -> FpNum<P> {
        if P % 4 == 3 {
            FpNum(Montgomery::<P>::from_u128(P - 1))
        } else if P % 8 == 3 || P % 8 == 5 {
            FpNum(Montgomery::<P>::from_u128(2))
        } else {
            let mut res = 0;
            let mut i = 0;
            while i < P {
                let a = standard_affine_shift(P, i);
                if intpow::<P>(a, (P - 1) / 2) == P - 1 {
                    res = a;
                    break;
                }
                i += 1;
            }
            FpNum(Montgomery::<P>::from_u128(res))
        }
    }

    /// Returns the underlying Montgomery number.
    pub fn montgomery(&self) -> Montgomery<P> {
        self.0
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
    const ONE: Self = FpNum(Montgomery::from_u128(1));
    const SIZE: u128 = P - 1;

    fn multiply(&self, other: &FpNum<P>) -> FpNum<P> {
        FpNum(self.0 * other.0)
    }
}

impl<const P: u128> From<u128> for FpNum<P> {
    fn from(value: u128) -> FpNum<P> {
        FpNum(Montgomery::<P>::from(value))
    }
}
impl<const P: u128> From<FpNum<P>> for u128 {
    fn from(value: FpNum<P>) -> u128 {
        u128::from(value.0)
    }
}
impl<const P: u128> From<&FpNum<P>> for u128 {
    fn from(value: &FpNum<P>) -> u128 {
        u128::from(value.0)
    }
}

impl<const P: u128> Neg for FpNum<P> {
    type Output = FpNum<P>;
    fn neg(self) -> FpNum<P> {
        FpNum(-self.0)
    }
}

impl<const P: u128> Add<Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, other: Self) -> FpNum<P> {
        FpNum(&self.0 + &other.0)
    }
}

impl<const P: u128> AddAssign<Self> for FpNum<P> {
    fn add_assign(&mut self, other: Self) {
        *self = FpNum(&self.0 + &other.0)
    }
}

impl<const P: u128> AddAssign<&Self> for FpNum<P> {
    fn add_assign(&mut self, other: &Self) {
        *self = FpNum(&self.0 + &other.0)
    }
}

impl<const P: u128> Sub<Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, other: Self) -> FpNum<P> {
        FpNum(&self.0 - &other.0)
    }
}

impl<const P: u128> SubAssign<&Self> for FpNum<P> {
    fn sub_assign(&mut self, other: &Self) {
        *self = FpNum(&self.0 - &other.0)
    }
}

impl<const P: u128> SubAssign<Self> for FpNum<P> {
    fn sub_assign(&mut self, other: Self) {
        *self = FpNum(&self.0 - &other.0)
    }
}

impl<const P: u128> Mul<&Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, other: &FpNum<P>) -> FpNum<P> {
        self.multiply(other)
    }
}

impl<const P: u128> Mul<Self> for FpNum<P> {
    type Output = FpNum<P>;

    fn mul(self, other: FpNum<P>) -> FpNum<P> {
        self.multiply(&other)
    }
}

impl<const P: u128> MulAssign<&Self> for FpNum<P> {
    fn mul_assign(&mut self, other: &Self) {
        *self = self.multiply(other)
    }
}

impl<const P: u128> MulAssign<Self> for FpNum<P> {
    fn mul_assign(&mut self, other: Self) {
        *self = self.multiply(&other);
    }
}

impl<const P: u128> Mul<FpNum<P>> for u128 {
    type Output = FpNum<P>;
    fn mul(self, other: FpNum<P>) -> FpNum<P> {
        FpNum::from(self).multiply(&other)
    }
}

impl<const P: u128> PartialEq<u128> for FpNum<P> {
    fn eq(&self, other: &u128) -> bool {
        u128::from(self) == *other
    }
}

impl<const P: u128> PartialEq<FpNum<P>> for u128 {
    fn eq(&self, other: &FpNum<P>) -> bool {
        *self == u128::from(other)
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
        assert_eq!(3, x);
        let five = FpNum::<7>::from(5);
        x = x.multiply(&five);
        assert_eq!(1, x);
        assert_eq!(5, five);
    }

    #[test]
    fn squares() {
        let mut x = FpNum::<7>::from(3);
        x *= x;
        assert_eq!(2, x);

        let mut x = FpNum::<7>::from(2);
        x *= x;
        assert_eq!(4, x);
    }

    #[test]
    fn powers_up() {
        let mut x = FpNum::<7>::from(2);
        x = x.pow(5);
        assert_eq!(4, x);

        let mut x = FpNum::<7>::from(3);
        x = x.pow(3);
        assert_eq!(6, x);

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
}
