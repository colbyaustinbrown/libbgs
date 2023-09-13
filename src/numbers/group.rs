use std::fmt;

use crate::numbers::Factorization;
use crate::numbers::SylowDecomposable;
use crate::util::*;

/// Types that represent the elements of a group.
/// In order for a type to represent the elements of the group, the type must satisfy these axioms:
/// * The type has a binary operator (`multiply`).
/// * The type must have a unique element given by `one()`.
/// * Event element of the type must have element given by `inverse()`.
pub trait GroupElem: Clone + PartialEq + Eq + fmt::Debug {
    /// Gets the size of the group this element belongs to.
    fn size() -> u128;
    
    /// Gets the unique identity element of this group.
    fn one() -> Self;

    /// True if this element is the multiplicative identity; false otherwise.
    fn is_one(&self) -> bool;

    /// Returns the product of two elements under the group binary operator.
    /// If you implement this trait, you must guarantee that the operation is associative; that is,
    /// `a.multiply(b.multiply(c, &g), &g) == a.multiply(b, &g).multiply(c, &g)`.
    fn multiply(&self, other: &Self) -> Self;

    /// Returns this element multiplied by itself.
    /// If you implement this trait, you must guarantee `x.square() == x.multiply(x)` for all `x`.
    fn square(&self) -> Self;

    /// Raises this element to the power of `n`.
    /// If you override this trait, you must guarantee that `x.pow(2) == x.square()` for all `x`.
    fn pow(&self, mut n: u128) -> Self {
        let mut y = Self::one();
        let mut res = self.clone();
        while n > 1 {
            // println!("{n} {:?} {:?}", &self, y);
            if n % 2 == 1 {
                y = y.multiply(&res);
            }
            res = res.square();
            n >>= 1;
        }
        res.multiply(&y)
    }

    /// Returns the multiplicative inverse of this element.
    /// If you implement this trait, you must guarantee `x.inverse().multiply(x)` and
    /// `x.multiply(x.inverse())` both evaluate to the unique identity element.
    fn inverse(&self) -> Self {
        self.pow(Self::size() - 1)
    }

    /// Returns the order of this element, that is, the smallest positive power `p` for which
    /// `a.pow(p, &g).is_one(&g)` returns True.
    fn order<const L: usize>(&self, parent_size: &Factorization<L>) -> Factorization<L> {
        let mut prime_powers = [(0, 0); L];
        for i in 0..L {
            let mut x = self.clone();
            for j in 0..L {
                if j == i {
                    continue;
                }
                x = x.pow(parent_size.factor(j));
            }

            let mut r = 0;
            while !x.is_one() {
                x = x.pow(parent_size[i].0);
                r += 1;
            }
            prime_powers[i] = (parent_size[i].0, r)
        }
        Factorization::new(prime_powers)
    }

    /// Returns the number of elements of a particular order.
    /// The argument is the powers of the prime factors of the group's order.
    fn count_elements_of_order<S, const L: usize>(ds: &[u128; L]) -> u128
    where
        // note: we only use the Factor trait here, but we require SylowDecomposable because this
        // count is only valid for finite cyclic groups.
        Self: SylowDecomposable<S, L>
    {
        let mut total = 1;
        for (d, (p, t)) in ds.iter().zip(Self::FACTORS.prime_powers()) {
            if *d > *t {
                return 0;
            } else if *d > 0 { 
                total *= intpow(*p, *d, 0) - intpow(*p, *d - 1, 0)
            }
        }
        total
    }
}
