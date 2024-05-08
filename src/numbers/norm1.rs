use std::ops::*;

use crate::numbers::*;

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

impl<const P: u128> Mul<Self> for Norm1<P> {
    type Output = Norm1<P>;
    fn mul(self, other: Self) -> Norm1<P> {
        Norm1(self.0 * other.0)
    }
}
