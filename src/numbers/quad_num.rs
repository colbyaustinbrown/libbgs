use std::ops::*;

use crate::numbers::*;

/// An element of the finite field $\mathbb{F}\_{p^2}$.
/// The number $x$ is represented as $x = a_0 + a_1\sqrt{r}$, where $r$ is
/// a basis element fixed at compile time.
/// See Lubeck, Frank. (2003). "Standard generators of finite fields and their cyclic subgroups."
/// Journal of Symbolic Computation (117) 51-67.
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
    const SIZE: u128 = P*P - 1;

    fn multiply(&self, other: &QuadNum<P>) -> QuadNum<P> {
        let a0 = self.0.multiply(&other.0) + self.1.multiply(&other.1).multiply(&QuadNum::<P>::R);
        let a1 = self.1.multiply(&other.0) + self.0.multiply(&other.1);

        QuadNum(a0, a1)
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
}
