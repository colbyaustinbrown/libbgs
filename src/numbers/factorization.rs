use std::ops::Index;

use crate::numbers::FactorStream;
use crate::util::intpow;

/// A prime power decomposition of a positive integer.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization<const L: usize> {
    value: u128,
    prime_powers: [(u128, u128); L],
}


/// Types that have a size or order which can be expressed as a product of prime powers.
/// The type parameter `S` is a phantom type to allow users of this library to provide their own
/// factorizations for `FpNum<P>`, `QuadNum<P>`, etc. for arbitrary `P`.
/// The type `S` must come from the crate implementing the factorization to satisfy Rust's
/// coherence and orphan rules, but does not affect the memory layout of any instance of an object
/// implementing `Factored`.
/// See the Rust documentation on [coherence and the orphan rule](https://doc.rust-lang.org/reference/items/implementations.html#trait-implementation-coherence).
pub trait Factored<S, const L: usize> {
    /// The prime factorization of this object.
    const FACTORS: Factorization<L>;
}

impl<const L: usize> Factorization<L> {
    /// Creates a new factorization from the given prime powers.
    pub const fn new(prime_powers: [(u128, u128); L]) -> Factorization<L> {
        let mut value = 1;
        let mut i = 0;
        while i < L {
            value *= intpow(prime_powers[i].0, prime_powers[i].1, 0);
            i += 1;
        }
        Factorization {
            value,
            prime_powers,
        }
    }

    /// Returns an `Iterator` yielding the maximal divisors of the number represented by this
    /// `Factorization` below `l`; that is, if this `Factorization` represents the number $n$, then
    /// this `Iterator` yields vectors representing all elements $d$ satisfying these properties:
    /// * $d | n$,
    /// * $d \leq l$, and
    /// * there are no elements $k | n$ such that $d | k$ and $k \leq l$.
    ///
    /// The iterator cannot outlive the `Factorization`, although the vectors yielded by it may.
    pub fn maximal_divisors(&self, l: u128) -> impl Iterator<Item = [u128; L]> + '_ {
        FactorStream::new(&self.prime_powers, l, true)
    }

    /// True if there the factorization represents 1.
    /// False otherwise.
    pub fn is_empty(&self) -> bool {
        self.prime_powers.is_empty()
    }

    /// Returns the prime power factor represented by prime number `i`, $p_i^{t_i}$.
    /// This method will `panic` if `i` is out of bounds.
    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    /// Gets the prime powers as an array.
    /// The first element of each entry is the prime, and the second is the power.
    pub fn prime_powers(&self) -> &[(u128, u128); L] {
        &self.prime_powers
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub fn value(&self) -> u128 {
        self.value
    }
}

impl<const L: usize> Index<usize> for Factorization<L> {
    type Output = (u128, u128);

    fn index(&self, index: usize) -> &(u128, u128) {
        &self.prime_powers[index]
    }
}
