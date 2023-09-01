use std::ops::*;

use crate::numbers::*;
use crate::util::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FpNum<const P: u128>(pub u128);

impl<const P: u128> FpNum<P> {
    pub fn int_sqrt(&self) -> Option<FpNum<P>> {
        if self.0 == 0 {
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
            let z = standard_affine_shift(P, i);
            if legendre(z, P) == (P - 1) {
                break z;
            }
            i += 1;
        };
        let mut c = FpNum::from(z).pow(q);
        let mut r = self.pow((q + 1) / 2);
        let mut t = self.pow(q);
        let mut m = s;

        loop {
            if t == 1 { return Some(r); }
            let mut temp = t;
            let mut i = 0;
            while temp != 1 {
                temp = temp.square();
                i += 1;
            }
            if i == m { return None; }
            let b = c.pow(1 << (m - i - 1));

            r *= b;
            t *= b.square();
            c = b.square();
            m = i;
        }
    }

    pub const fn find_nonresidue(p: u128) -> u128 {
        if p % 4 == 3 {
            p - 1
        } else if p % 8 == 3 || p % 8 == 5 {
            2
        } else {
            let mut res = 0;
            let mut i = 0;
            while i < p {
                let a = standard_affine_shift(p, i);
                if intpow(a, (p - 1) / 2, p) == p - 1 {
                    res = a;
                    break;
                }
                i += 1;
            }
            res
        }
    }
}

impl<S, const P: u128, const L: usize> SylowDecomposable<S, L> for FpNum<P> 
where
    FpNum<P>: Factored<S, L>
{
    type Elem = FpNum<P>;
    fn find_sylow_generator(i: usize) -> FpNum<P> {
        match <Self as Factored<S, L>>::FACTORS[i] {
            (2, 1) => FpNum::from(FpNum::<P>::size()),
            (p, t) => (1..FpNum::<P>::size())
                .map(|j| FpNum::from(standard_affine_shift(P, j)))
                .find_map(|c| FpNum::is_sylow_generator(&c, (p, t)))
                .unwrap(),
        }
    }
}

impl<const P: u128> GroupElem for FpNum<P> {
    fn is_one(&self) -> bool {
        self.0 == 1
    }

    fn multiply(&self, other: &FpNum<P>) -> FpNum<P> {
        *self * *other
    }

    fn square(&self) -> FpNum<P> {
        *self * *self
    }

    fn size() -> u128 {
        P - 1
    }

    fn one() -> FpNum<P> {
        FpNum::from(1)
    }
}

impl<const P: u128> From<u128> for FpNum<P> {
    fn from(value: u128) -> FpNum<P> {
        FpNum(value)
    }
}
impl<const P: u128> From<FpNum<P>> for u128 {
    fn from(value: FpNum<P>) -> u128 {
        value.0
    }
}
impl<const P: u128> From<&FpNum<P>> for u128 {
    fn from(value: &FpNum<P>) -> u128 {
        value.0
    }
}

impl<const P: u128> Add<Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn add(self, other: Self) -> FpNum<P> {
        FpNum((self.0 + other.0) % P)
    }
}

impl<const P: u128> AddAssign<Self> for FpNum<P> {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
        self.0 %= P;
    }
}

impl<const P: u128> AddAssign<&Self> for FpNum<P> {
    fn add_assign(&mut self, other: &Self) {
        self.0 += other.0;
        self.0 %= P;
    }
}

impl<const P: u128> AddAssign<u128> for FpNum<P> {
    fn add_assign(&mut self, other: u128) {
        self.0 += other;
        self.0 %= P;
    }
}

impl<const P: u128> Sub<Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn sub(self, other: Self) -> FpNum<P> {
        let mut value = self.0;
        while value < other.0 {
            value += P;
        }
        value -= other.0;
        FpNum(value)
    }
}

impl<const P: u128> SubAssign<&Self> for FpNum<P> {
    fn sub_assign(&mut self, other: &Self) {
        while self.0 < other.0 {
            self.0 += P;
        }
        self.0 -= other.0;
    }
}

impl<const P: u128> SubAssign<Self> for FpNum<P> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<const P: u128> SubAssign<u128> for FpNum<P> {
    fn sub_assign(&mut self, other: u128) {
        while self.0 < other {
            self.0 += P;
        }
        self.0 -= other;
    }
}

impl<const P: u128> Mul<&Self> for FpNum<P> {
    type Output = FpNum<P>;
    fn mul(self, other: &FpNum<P>) -> FpNum<P> {
        FpNum(long_multiply(self.0, other.0, P))
    }
}

impl<const P: u128> Mul<Self> for FpNum<P> {
    type Output = FpNum<P>;

    fn mul(self, other: FpNum<P>) -> FpNum<P> {
        FpNum(long_multiply(self.0, other.0, P))
    }
}

impl<const P: u128> MulAssign<&Self> for FpNum<P> {
    fn mul_assign(&mut self, other: &Self) {
        self.0 = long_multiply(self.0, other.0, P);
    }
}

impl<const P: u128> MulAssign<Self> for FpNum<P> {
    fn mul_assign(&mut self, other: Self) {
        self.0 = long_multiply(self.0, other.0, P);
    }
}

impl<const P: u128> Mul<FpNum<P>> for u128 {
    type Output = FpNum<P>;
    fn mul(self, other: FpNum<P>) -> FpNum<P> {
        FpNum(long_multiply(self, other.0, P))
    }
}

impl<const P: u128> PartialEq<u128> for FpNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.0 == *other
    }
}

impl<const P: u128> std::fmt::Display for FpNum<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        self.0.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    /*
    #[test]
    fn one_is_one() {
        let p = FpStar::<7> {};
        let one = p.one();
        assert!(one.is_one(&p));
    }

    #[test]
    fn multiplies() {
        let p = FpStar::<7> {};
        let mut x = FpNum::from(3);
        assert_eq!(3, x.0);
        let five = FpNum::from(5);
        x = x.multiply(&five, &p);
        assert_eq!(1, x.0);
        assert_eq!(5, five.0);
    }

    #[test]
    fn squares() {
        let mut x = FpNum::<7>::from(3);
        x *= x;
        assert_eq!(2, x.0);

        let mut x = FpNum::<7>::from(2);
        x *= x;
        assert_eq!(4, x.0);
    }

    #[test]
    fn powers_up() {
        let p = FpStar::<7> {};
        let mut x = FpNum::from(2);
        x = x.pow(5, &p);
        assert_eq!(4, x.0);

        let mut x = FpNum::from(3);
        x = x.pow(3, &p);
        assert_eq!(6, x.0);

        let mut x = FpNum::from(5);
        x = x.pow(6, &p);
        assert!(x.is_one(&p));
    }

    #[test]
    fn powers_up_big() {
        let p = FpStar::<BIG_P> {};
        let mut x = FpNum::from(3);
        x = x.pow(BIG_P - 1, &FpStar::<BIG_P> {});
        println!("x is {x:?}");
        assert!(x.is_one(&p));
    }

    #[test]
    fn sylow_one_is_one() {
        let p = Factorization::new(vec![(2, 2), (3, 1)]);
        let g = SylowDecomp::new(&FpStar::<13> {}, p.clone());
        let one = g.one();
        assert!(one.is_one(&g));
    }

    #[test]
    fn sylow_finds_generators() {
        let p = Factorization::new(vec![(2, 2), (7, 1)]);
        let g = SylowDecomp::new(&FpStar::<29> {}, p.clone());
        for i in 0..g.generators().len() {
            let gen = &g.generators()[i];
            let d = g.factors().factor(i);
            test_is_generator_small::<FpStar<29>>(gen, d, &FpStar::<29> {});
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fact = Factorization::new(vec![
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
        let g = SylowDecomp::new(&FpStar::<BIG_P> {}, fact.clone());
        for i in 0..g.generators().len() {
            let gen = &g.generators()[i];
            let d = fact[i];
            test_is_generator_big::<FpStar<BIG_P>>(gen, d, &FpStar::<BIG_P> {});
        }
    }

    #[test]
    fn sylow_order() {
        let fact = Factorization::new(vec![(2, 2), (3, 1)]);
        let g = SylowDecomp::new(&FpStar::<13> {}, fact.clone());
        for i in 1..13 {
            let mut x = SylowElem::new(vec![i % 4, i % 3]);
            x = x.pow(x.order(&g).value(), &g);
            assert!(x.is_one(&g));
        }
    }

    #[test]
    fn sylow_order_big() {
        let fact = Factorization::new(vec![
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
        let g = SylowDecomp::new(&FpStar::<BIG_P> {}, fact);
        let n = 123456789;
        let mut x = SylowElem::new(
            g.factors()
                .as_array()
                .iter()
                .map(|(p, d)| n % intpow(*p, *d, 0))
                .collect(),
        );
        let or = x.order(&g).value();
        x = x.pow(or, &g);
        assert!(x.is_one(&g));
    }

    #[test]
    fn calculates_square_roots() {
        let p = FpStar::<13> {};
        let mut nonresidues = 0;
        for x in (1..13).map(|i| FpNum::from(i)) {
            match x.int_sqrt() {
                None => {
                    nonresidues += 1;
                }
                Some(mut y) => {
                    y = y.pow(2, &p);
                    assert_eq!(x, y);
                }
            }
        }
        assert_eq!(nonresidues, 6);
    }

    #[test]
    fn inverses() {
        let p = FpStar::<13> {};
        for i in 2..13 {
            let mut x = FpNum::from(i);
            let y = x.clone();
            println!("{x:?}");
            x = x.inverse(&p);
            println!("{x:?}");
            assert!(!x.is_one(&p));
            x = x.multiply(&y, &p);
            assert!(x.is_one(&p));
        }
    }
    */
}
