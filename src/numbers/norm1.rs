use std::ops::*;

use crate::numbers::*;
use crate::util::*;

/// An element $x \in \mathbb{F}\_{p^2}$ with norm $N\_{\mathbb{F}\_{p^2}/\mathbb{F}\_p}(x) = 1$.
/// This is equivalent to $x$ having order dividing $p + 1$.
#[derive(PartialEq, Eq, Clone, Copy)]
pub struct Norm1<const P: u128>(QuadNum<P>);

impl<const P: u128> GroupElem for Norm1<P> {
    const ONE: Self = Norm1(QuadNum::ONE);
    const SIZE: u128 = P + 1;

    fn multiply(&self, other: &Norm1<P>) -> Norm1<P> {
        *self * *other
    }
}

impl<S, const P: u128> SylowDecomposable<S> for Norm1<P>
where
    Norm1<P>: FactoredSize<S>,
{
    fn find_sylow_generator(i: usize) -> Norm1<P> {
        (1..P * 2)
            .map(|i| {
                let j = standard_affine_shift(P * 2, i);
                let p = QuadNum::steinitz(j);
                p.pow(P - 1)
            })
            .filter(|c| *c != QuadNum::ZERO)
            .find_map(|c| Norm1::is_sylow_generator(&Norm1(c), Self::FACTORS[i]))
            .unwrap()
    }
}

impl<const P: u128> From<Norm1<P>> for QuadNum<P> {
    fn from(src: Norm1<P>) -> QuadNum<P> {
        src.0 
    }
}

impl<const P: u128> Add<Self> for Norm1<P> {
    type Output = QuadNum<P>;
    fn add(self, other: Self) -> QuadNum<P> {
        self.0 + other.0
    }
}

impl<const P: u128> Sub<Self> for Norm1<P> {
    type Output = QuadNum<P>;
    fn sub(self, other: Self) -> QuadNum<P> {
        self.0 - other.0
    }
}

impl<const P: u128> Mul<Self> for Norm1<P> {
    type Output = Norm1<P>;
    fn mul(self, other: Self) -> Norm1<P> {
        Norm1(self.0 * other.0)
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

    impl FactoredSize<Phantom> for Norm1<7> {
        const FACTORS: Factorization = Factorization::new(&[(2, 3)]);
    }

    impl FactoredSize<Phantom> for Norm1<17> {
        const FACTORS: Factorization = Factorization::new(&[(2, 1), (3, 2)]);
    }

    impl_factors!(Phantom, 41);

    #[test]
    fn powers_up_big() {
        let mut x = QuadNum::<BIG_P>(FpNum::from(3), FpNum::from(5));
        x = x.pow(BIG_P - 1);
        x = x.pow(BIG_P + 1);
        assert!(x == QuadNum::ONE);
    }

    #[test]
    fn sylow_finds_generators() {
        let g = SylowDecomp::<Phantom, 2, Norm1<17>>::new();
        for i in 0..2 {
            let gen = &g.generator(i);
            let d = SylowElem::<Phantom, 2, Norm1<17>>::FACTORS.factor(i);
            test_is_generator_small::<Phantom, 2, Norm1<17>>(*gen, d as usize);
        }
    }

    #[test]
    fn sylow_finds_generators_2() {
        let g = SylowDecomp::<Phantom, 3, Norm1<41>>::new();
        for i in 0..3 {
            let gen = g.generator(i);
            assert!(gen.0 != QuadNum(FpNum::from(0), FpNum::from(0)));
            let d = SylowElem::<Phantom, 3, Norm1<41>>::FACTORS.factor(i);
            test_is_generator_small::<Phantom, 2, Norm1<41>>(gen, d as usize);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let g = SylowDecomp::<Phantom, 11, Norm1<BIG_P>>::new();
        for i in 0..11 {
            let gen = g.generator(i);
            let d = SylowElem::<Phantom, 11, Norm1<BIG_P>>::FACTORS[i];
            test_is_generator_big::<Phantom, 11, Norm1<BIG_P>>(gen, d);
        }
    }
}
