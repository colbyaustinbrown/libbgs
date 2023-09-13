use std::ops::*;

use either::*;

use crate::numbers::*;
use crate::util::*;

/// The finite field of size `P^2`. Isomorphic to $\mathbb{Z} / P^2\mathbb{Z}$.
/// Each `QuadField` has a fixed quadratic nonresidue `r` used as a basis element for the numbers
/// outside of the prime subfield.
#[derive(PartialEq, Eq, Debug)]
pub struct QuadField<const P: u128> {
    r: u128,
}

/// An integer modulo `P^2`. An element $x$ is represented as $x = a_0 + a_1\sqrt{r}$, where $r$ is
/// the fixed basis element.
///
/// The association between a `QuadNum` instance and a `QuadField` (hence, the fixed basis element
/// `r`), is implicit.
/// You must ensure that operations on `QuadNum` instances are associated to the same choice of
/// `r`.
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub struct QuadNum<const P: u128>(
    /// The value $a_0$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub u128, 
    /// The value $a_1$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub u128,
);

impl<S, const P: u128, const L: usize> SylowDecomposable<S, L> for QuadNum<P>
where
    QuadNum<P>: Factor<S, L>,
{
    fn find_sylow_generator(i: usize) -> QuadNum<P> {
        let pow = P - 1;
        // should be self.p * self.p, but maybe this works?
        (1..P * 2)
            .map(|i| {
                let j = standard_affine_shift(P * 2, i);
                let p = QuadNum::steinitz(j);
                p.pow(pow)
            })
            .find_map(|c| QuadNum::is_sylow_generator(&c, <Self as Factor<S, L>>::FACTORS[i]))
            .unwrap()
    }
}

impl<const P: u128> QuadNum<P> {
    /// True if this number is zero; false otherwise.
    pub fn is_zero(&self) -> bool {
        self.0 == 0 && self.1 == 0
    }

    /// Returns the Steinitz element of $\mathbb{F}\_{p^2}$ with index `i`.
    pub fn steinitz(i: u128) -> QuadNum<P> {
        QuadNum::from((i % P, i / P))
    }

    /// Calculates the square root of an integer modulo `P`, casting to an `FpNum<P>` if `x` is a
    /// quadratic residue.
    /// Returns a `Left` `QuadNum<P>` if `x` is a quadratic nonresidue, or a `Right` `FpNum<P>` if
    /// `x` is a quadratic residue (including 0).
    pub fn int_sqrt_either(x: u128) -> Either<QuadNum<P>, FpNum<P>> {
        let mut x = FpNum::from(x);
        if let Some(y) = x.int_sqrt() {
            return Right(y);
        }

        let r = FpNum::from(Self::R).inverse();
        x = x.multiply(&r);
        let a1 = x.int_sqrt().unwrap();
        Left(QuadNum(0, a1.into()))
    }

    /// Calculates the square root af in integer modulo `P`.
    pub fn int_sqrt(x: u128) -> QuadNum<P> {
        Self::int_sqrt_either(x).left_or_else(|n| QuadNum::from((n.into(), 0)))
    }

    /// The basis element for the numbers outside of the prime subfield.
    pub const R: u128 = FpNum::<P>::find_nonresidue(P);
}

impl<const P: u128> GroupElem for QuadNum<P> {
    fn is_one(&self) -> bool {
        self.0 == 1 && self.1 == 0
    }

    fn multiply(&self, other: &QuadNum<P>) -> QuadNum<P> {
        let mut a0 = long_multiply(self.0, other.0, P) + long_multiply(self.1, long_multiply(other.1, QuadNum::<P>::R, P), P);
        let mut a1 = long_multiply(self.1, other.0, P) + long_multiply(self.0, other.1, P);
        if a0 >= P {
            a0 -= P;
        }
        if a1 >= P {
            a1 -= P;
        }
        QuadNum(a0, a1)
    }

    fn size() -> u128 {
        P + 1
    }

    fn one() -> QuadNum<P> {
        QuadNum(1, 0)
    }
}

impl<const P: u128> PartialEq<u128> for QuadNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.0 == *other && self.1 == 0
    }
}

impl<const P: u128> From<(u128, u128)> for QuadNum<P> {
    fn from(value: (u128, u128)) -> QuadNum<P> {
        QuadNum(value.0, value.1)
    }
}

impl<const P: u128> Add<Self> for QuadNum<P> {
    type Output = QuadNum<P>;
    fn add(self, other: Self) -> QuadNum<P> {
        let mut a0 = self.0 + other.0;
        let mut a1 = self.1 + other.1;
        if a0 >= P {
            a0 -= P;
        }
        if a1 >= P {
            a1 -= P;
        }
        QuadNum(a0, a1)
    }
}

impl<const P: u128> AddAssign<Self> for QuadNum<P> {
    fn add_assign(&mut self, other: Self) {
        self.0 = self.0 + other.0;
        self.1 = self.1 + other.1;
        if self.0 >= P {
            self.0 -= P;
        }
        if self.1 >= P {
            self.1 -= P;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    impl Factor<Phantom, 1> for QuadNum<7> {
        const FACTORS: Factorization<1> = Factorization::new([(2, 3)]);
    }

    impl Factor<Phantom, 2> for QuadNum<17> {
        const FACTORS: Factorization<2> = Factorization::new([(2, 1), (3, 2)]);
    }

    impl Factor<Phantom, 11> for QuadNum<BIG_P> {
        const FACTORS: Factorization<11> = Factorization::new([
            (2, 4),
            (3, 1),
            (5, 2),
            (11, 2),
            (17, 1),
            (19, 1),
            (23, 1),
            (97, 1),
            (757, 1),
            (1453, 1),
            (8689, 1),
        ]);
    }

    #[test]
    fn one_is_one() {
        let one = QuadNum::<7>::one();
        assert!(one.is_one());
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        for i in 2..7 {
            assert_ne!((i * i) % 7, QuadNum::<7>::R);
        }
    }

    #[test]
    fn powers_up() {
        let mut x = QuadNum::<7>::from((3, 4));
        x = x.pow(48);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let mut x = QuadNum::<BIG_P>(3, 5);
        x = x.pow(BIG_P - 1);
        x = x.pow(BIG_P + 1);
        assert!(x.is_one());
    }

    #[test]
    fn finds_sqrt() {
        for i in 3..1003 {
            let mut x = QuadNum::<BIG_P>::int_sqrt(i);
            assert_ne!(x, i);
            x = x.multiply(&x);
            assert_eq!(x, i);
        }
    }

    #[test]
    fn sylow_finds_generators() {
        let g = SylowDecomp::<Phantom, 2, QuadNum<17>>::new();
        for i in 0..g.generators().len() {
            let gen = &g.generators()[i];
            let d = SylowElem::<Phantom, 2, QuadNum<17>>::FACTORS.factor(i);
            test_is_generator_small(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let g = SylowDecomp::<Phantom, 11, QuadNum<BIG_P>>::new();
        for i in 0..g.generators().len() {
            let gen = &g.generators()[i];
            let d = SylowElem::<Phantom, 11, QuadNum<BIG_P>>::FACTORS.prime_powers()[i];
            test_is_generator_big(gen, d);
        }
    }
}
