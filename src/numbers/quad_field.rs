use std::ops::*;

use either::*;

use crate::numbers::*;
use libbgs_util::*;

/// An integer modulo `P^2`. An element $x$ is represented as $x = a_0 + a_1\sqrt{r}$, where $r$ is
/// the fixed basis element.
/// See Lubeck, Frank. (2003). "Standard generators of finite fields and their cyclic subgroups."
/// Journal of Symbolic Computation (117) 51-67.
#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash)]
pub struct QuadNum<const P: u128>(
    /// The value $a_0$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub FpNum<P>,
    /// The value $a_1$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub FpNum<P>,
);

impl<const P: u128> QuadNum<P> {
    /// The basis element for the numbers outside of the prime subfield.
    pub const R: FpNum<P> = FpNum::<P>::find_nonresidue();

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
    pub fn int_sqrt_either(mut x: FpNum<P>) -> Either<QuadNum<P>, FpNum<P>> {
        if let Some(y) = x.int_sqrt() {
            return Right(y);
        }

        let r = Self::R.inverse();
        x = x.multiply(&r);
        let a1 = x.int_sqrt().unwrap();
        Left(QuadNum(FpNum::from(0), a1))
    }

    /// Calculates the square root af in integer modulo `P`.
    pub fn int_sqrt(x: FpNum<P>) -> QuadNum<P> {
        Self::int_sqrt_either(x).left_or_else(|n| QuadNum::from((n.into(), 0)))
    }
}

impl<const P: u128> GroupElem for QuadNum<P> {
    const ONE: Self = QuadNum(
        FpNum(Montgomery::from_u128(1)),
        FpNum(Montgomery::from_u128(0)),
    );
    const SIZE: u128 = P + 1;

    fn multiply(&self, other: &QuadNum<P>) -> QuadNum<P> {
        let a0 = self.0.multiply(&other.0) + self.1.multiply(&other.1).multiply(&QuadNum::<P>::R);
        let a1 = self.1.multiply(&other.0) + self.0.multiply(&other.1);

        QuadNum(a0, a1)
    }
}

impl<S, const P: u128> SylowDecomposable<S> for QuadNum<P>
where
    QuadNum<P>: Factor<S>,
{
    fn find_sylow_generator(i: usize) -> QuadNum<P> {
        (1..P * 2)
            .map(|i| {
                let j = standard_affine_shift(P * 2, i);
                let p = QuadNum::steinitz(j);
                p.pow(P - 1)
            })
            .find_map(|c| {
                QuadNum::is_sylow_generator(&c, Self::FACTORS[i])
            })
            .unwrap()
    }
}

impl<const P: u128> PartialEq<u128> for QuadNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.0 == *other && self.1 == 0
    }
}

impl<const P: u128> From<FpNum<P>> for QuadNum<P> {
    fn from(value: FpNum<P>) -> QuadNum<P> {
        QuadNum(value, FpNum::from(0))
    }
}

impl<const P: u128> From<(u128, u128)> for QuadNum<P> {
    fn from(value: (u128, u128)) -> QuadNum<P> {
        QuadNum(FpNum::from(value.0), FpNum::from(value.1))
    }
}

impl<const P: u128> Add<Self> for QuadNum<P> {
    type Output = QuadNum<P>;
    fn add(self, other: Self) -> QuadNum<P> {
        let a0 = self.0 + other.0;
        let a1 = self.1 + other.1;
        QuadNum(a0, a1)
    }
}

impl<const P: u128> Sub<Self> for QuadNum<P> {
    type Output = QuadNum<P>;
    fn sub(self, other: Self) -> QuadNum<P> {
        let a0 = self.0 - other.0;
        let a1 = self.1 - other.1;
        QuadNum(a0, a1)
    }
}

impl<const P: u128> AddAssign<Self> for QuadNum<P> {
    fn add_assign(&mut self, other: Self) {
        self.0 = self.0 + other.0;
        self.1 = self.1 + other.1;
    }
}

impl<const P: u128> Mul<Self> for QuadNum<P> {
    type Output = QuadNum<P>;
    fn mul(self, other: Self) -> QuadNum<P> {
        self.multiply(&other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    impl_factors!(Phantom, 1_000_000_000_000_000_124_399);

    impl Factor<Phantom> for QuadNum<7> {
        const FACTORS: Factorization = Factorization::new(&[(2, 3)]);
    }

    impl Factor<Phantom> for QuadNum<17> {
        const FACTORS: Factorization = Factorization::new(&[(2, 1), (3, 2)]);
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        for i in 2..7 {
            assert_ne!((i * i) % 7, u128::from(QuadNum::<7>::R));
        }
    }

    #[test]
    fn powers_up() {
        let mut x = QuadNum::<7>::from((3, 4));
        x = x.pow(48);
        assert!(x == QuadNum::ONE);
    }

    #[test]
    fn powers_up_big() {
        let mut x = QuadNum::<BIG_P>::from((3, 5));
        x = x.pow(BIG_P - 1);
        x = x.pow(BIG_P + 1);
        assert!(x == QuadNum::ONE);
    }

    #[test]
    fn finds_sqrt() {
        for i in 3..1003 {
            let mut x = QuadNum::<BIG_P>::int_sqrt(FpNum::from(i));
            assert_ne!(x, i);
            x = x.multiply(&x);
            assert_eq!(x, i);
        }
    }

    #[test]
    fn sylow_finds_generators() {
        let g = SylowDecomp::<Phantom, 2, QuadNum<17>>::new();
        for i in 0..2 {
            let gen = &g.generator(i);
            let d = SylowElem::<Phantom, 2, QuadNum<17>>::FACTORS.factor(i);
            test_is_generator_small::<Phantom, 2, QuadNum<17>>(*gen, d as usize);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let g = SylowDecomp::<Phantom, 11, QuadNum<BIG_P>>::new();
        for i in 0..11 {
            let gen = g.generator(i);
            let d = SylowElem::<Phantom, 11, QuadNum<BIG_P>>::FACTORS[i];
            test_is_generator_big::<Phantom, 11, QuadNum<BIG_P>>(gen, d);
        }
    }
}
