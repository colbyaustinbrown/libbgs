use std::ops::*;

use crate::numbers::*;
use libbgs_util::*;

/// An element of the finite field $\mathbb{F}\_{p^2}$.
/// The number $x$ is represented as $x = a_0 + a_1\sqrt{r}$, where $r$ is
/// a basis element fixed at compile time.
/// See Lubeck, Frank. (2003). "Standard generators of finite fields and their cyclic subgroups."
/// Journal of Symbolic Computation (117) 51-67.
/// Note that the `SylowDecomposable` implementation for a `QuadNum` returns the decomposition for
/// the subgroup with $p + 1$ elements, not the full group $\mathbb{F}_{p^2}^\times$.
/// Also, `<QuadNum<P> as GroupElem>::SIZE == P + 1`, again refering to the subgroup.
/// For these reasons, this API is likely to change in the future to bring the definitions of `QuadNum<P> as
/// GroupElem` and the `SylowDecomp` instance in line with describing the full group.
#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash)]
pub struct QuadNum<const P: u128>(
    /// The value $a_0$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub(crate) FpNum<P>,
    /// The value $a_1$, when writing this `QuadNum` as $a_0 + a_1\sqrt{r}$.
    pub(crate) FpNum<P>,
);

impl<const P: u128> QuadNum<P> {
    /// The basis element for the numbers outside of the prime subfield.
    pub const R: FpNum<P> = FpNum::<P>::find_nonresidue();

    /// The constant zero.
    pub const ZERO: QuadNum<P> = QuadNum(FpNum::from_u128(0), FpNum::from_u128(0));

    /// True if this number is zero; false otherwise.
    pub const fn is_zero(&self) -> bool {
        self.0.raw() == FpNum::<P>::ZERO.raw() && self.1.raw() == FpNum::<P>::ZERO.raw()
    }

    /// Returns the Steinitz element of $\mathbb{F}\_{p^2}$ with index `i`.
    pub fn steinitz(i: u128) -> QuadNum<P> {
        QuadNum(FpNum::from(i % P), FpNum::from(i / P))
    }

    /// Returns the Steinitz element of $\mathbb{F}\_{p^2}$ with index `i`.
    /// This method is more expensive than `steinitz`, but can be used in `const` contexts.
    pub const fn steinitz_const(i: u128) -> QuadNum<P> {
        QuadNum(FpNum::from_u128(i % P), FpNum::from_u128(i / P))
    }
}

impl<const P: u128> GroupElem for QuadNum<P> {
    const ONE: Self = QuadNum(
        FpNum::ONE,
        FpNum::ZERO,
    );
    const SIZE: u128 = P + 1;

    fn multiply(&self, other: &QuadNum<P>) -> QuadNum<P> {
        let a0 = self.0.multiply(&other.0) + self.1.multiply(&other.1).multiply(&QuadNum::<P>::R);
        let a1 = self.1.multiply(&other.0) + self.0.multiply(&other.1);

        QuadNum(a0, a1)
    }

    fn inverse(&self) -> QuadNum<P> {
        if *self == QuadNum::ZERO {
            panic!("Attempted to take the multiplicative inverse of zero."); 
        }
        self.pow(P * P - 2)
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
            .filter(|c| *c != QuadNum::ZERO)
            .find_map(|c| QuadNum::is_sylow_generator(&c, Self::FACTORS[i]))
            .unwrap()
    }
}

impl<const P: u128> PartialEq<u128> for QuadNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.0 == FpNum::from(*other) && self.1 == FpNum::ZERO 
    }
}

impl<const P: u128> From<FpNum<P>> for QuadNum<P> {
    fn from(value: FpNum<P>) -> QuadNum<P> {
        QuadNum(value, FpNum::from(0))
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

    impl_factors!(Phantom, 41);

    #[test]
    fn calculates_r_as_nonresidue() {
        for i in 2..7 {
            assert_ne!((i * i) % 7, u128::from(QuadNum::<7>::R));
        }
    }

    #[test]
    fn powers_up() {
        let mut x = QuadNum::<7>(FpNum::from(3), FpNum::from(4));
        x = x.pow(48);
        assert!(x == QuadNum::ONE);
    }

    #[test]
    fn powers_up_big() {
        let mut x = QuadNum::<BIG_P>(FpNum::from(3), FpNum::from(5));
        x = x.pow(BIG_P - 1);
        x = x.pow(BIG_P + 1);
        assert!(x == QuadNum::ONE);
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
    fn sylow_finds_generators_2() {
        let g = SylowDecomp::<Phantom, 3, QuadNum<41>>::new();
        for i in 0..3 {
            let gen = g.generator(i);
            assert!(*gen != QuadNum(FpNum::from(0), FpNum::from(0)));
            let d = SylowElem::<Phantom, 3, QuadNum<41>>::FACTORS.factor(i);
            test_is_generator_small::<Phantom, 2, QuadNum<41>>(gen, d as usize);
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
