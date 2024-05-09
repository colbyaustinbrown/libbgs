use crate::numbers::FactoredSize;
use libbgs_util::*;

/// Types that represent the elements of a group.
/// In order for a type to represent the elements of the group, the type must satisfy these axioms:
/// * The type has a binary operator (`multiply`).
/// * The type must have a `const ONE`
/// * Every element of the type must have element given by `inverse()`.
pub trait GroupElem: Clone + Eq {
    /// The unique identity element of this group.
    const ONE: Self;

    /// 256 copies of `Self::ONE` in an array.
    #[deprecated(note = "To be replaced by inline `const` expressions once stabilized.")]
    const ONE_256: [Self; 256] = [Self::ONE; 256];

    /// The size of the group this element belongs to.
    const SIZE: u128;

    /// Returns the product of two elements under the group binary operator.
    /// If you implement this trait, you must guarantee that the operation is associative; that is,
    /// `a.multiply(b.multiply(c)) == a.multiply(b).multiply(c)`.
    fn multiply(&self, other: &Self) -> Self;

    /// Raises this element to the power of `n`.
    fn pow(&self, mut n: u128) -> Self {
        if n == 0 {
            return Self::ONE;
        }
        let mut y = Self::ONE;
        let mut res = self.clone();
        while n > 1 {
            if n % 2 == 1 {
                y = y.multiply(&res);
            }
            res = res.multiply(&res);
            n /= 2;
        }
        res.multiply(&y)
    }

    /// Returns the multiplicative inverse of this element.
    /// If you implement this trait, you must guarantee `x.inverse().multiply(x)` and
    /// `x.multiply(x.inverse())` both evaluate to `ONE`.
    fn inverse(&self) -> Self {
        self.pow(Self::SIZE - 1)
    }

    /// Returns the order of this element, that is, the smallest positive power `p` for which
    /// `a.pow(p).is_one()` returns True.
    fn order<S>(&self) -> u128
    where
        Self: FactoredSize<S>,
    {
        let mut res = 1;
        for i in 0..Self::FACTORS.len() {
            let mut x = self.clone();
            for j in 0..Self::FACTORS.len() {
                if j == i {
                    continue;
                }
                x = x.pow(Self::FACTORS.factor(j));
            }

            let mut r = 0;
            while x != Self::ONE {
                x = x.pow(Self::FACTORS[i].0);
                r += 1;
            }
            res *= intpow::<0>(Self::FACTORS[i].0, r);
        }
        res
    }
}
