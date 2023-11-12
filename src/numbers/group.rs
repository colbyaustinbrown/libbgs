use std::fmt;

use crate::numbers::SylowDecomposable;
use crate::numbers::{Factor, Length};
use libbgs_util::*;

/// Types that represent the elements of a group.
/// In order for a type to represent the elements of the group, the type must satisfy these axioms:
/// * The type has a binary operator (`multiply`).
/// * The type must have a unique element given by `one()`.
/// * Event element of the type must have element given by `inverse()`.
pub trait GroupElem: Clone + PartialEq + Eq + fmt::Debug {
    /// Gets the unique identity element of this group.
    const ONE: Self;

    ///! 256 copies of `Self::ONE` in an array.
    #[deprecated(note = "To be replaced by inline `const` expressions once stabilized.")]
    const ONE_256: [Self; 256] = [Self::ONE; 256];

    /// Gets the size of the group this element belongs to.
    const SIZE: u128;

    /// Returns the product of two elements under the group binary operator.
    /// If you implement this trait, you must guarantee that the operation is associative; that is,
    /// `a.multiply(b.multiply(c, &g), &g) == a.multiply(b, &g).multiply(c, &g)`.
    fn multiply(&self, other: &Self) -> Self;

    /// Raises this element to the power of `n`.
    fn pow(&self, mut n: u128) -> Self {
        let mut y = Self::ONE;
        let mut res = self.clone();
        while n > 1 {
            if n & 1 == 1 {
                y = y.multiply(&res);
            }
            res = res.multiply(&res);
            n >>= 1;
        }
        res.multiply(&y)
    }

    /// Returns the multiplicative inverse of this element.
    /// If you implement this trait, you must guarantee `x.inverse().multiply(x)` and
    /// `x.multiply(x.inverse())` both evaluate to the unique identity element.
    fn inverse(&self) -> Self {
        self.pow(Self::SIZE - 1)
    }

    /// Returns the order of this element, that is, the smallest positive power `p` for which
    /// `a.pow(p, &g).is_one(&g)` returns True.
    fn order<S>(&self) -> u128
    where
        Self: Factor<S>,
    {
        let mut res = 1;
        for i in 0..Self::LENGTH {
            let mut x = self.clone();
            for j in 0..Self::LENGTH {
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

    /// Returns the number of elements of a particular order.
    /// The argument is the powers of the prime factors of the group's order.
    fn count_elements_of_order<S, const L: usize>(ds: &[usize; L]) -> u128
    where
        // note: we only use the Factor trait here, but we require SylowDecomposable because this
        // count is only valid for finite cyclic groups.
        Self: SylowDecomposable<S>,
    {
        let mut total = 1;
        for (d, (p, t)) in ds.iter().zip(Self::FACTORS.prime_powers()) {
            if *d > *t {
                return 0;
            } else if *d > 0 {
                let tmp = intpow::<0>(*p, (*d - 1) as u128);
                total *= tmp * *p - tmp;
            }
        }
        total
    }
}
