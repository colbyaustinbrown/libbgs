use std::ops::Index;

use crate::numbers::FactorStream;
use crate::util::intpow;

/// A prime power decomposition of a positive integer.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    value: u128,
    prime_powers: &'static [(u128, usize)],
}


/// Types that have a size or order which can be expressed as a product of prime powers.
/// The type parameter `S` is a phantom type to allow users of this library to provide their own
/// factorizations for `FpNum<P>`, `QuadNum<P>`, etc. for arbitrary `P`.
/// The type `S` must come from the crate implementing the factorization to satisfy Rust's
/// coherence and orphan rules, but does not affect the memory layout of any instance of an object
/// implementing `Factor`.
/// See the Rust documentation on [coherence and the orphan rule](https://doc.rust-lang.org/reference/items/implementations.html#trait-implementation-coherence).
///
/// This type can only hold factors up to `2^126`. Behavior may be undefined for factorizations
/// which contain larger factors.
pub trait Factor<S> {
    /// The number of prime factors in the factorization.
    const LEN: usize = Self::FACTORS.len();
    /// The prime factorization of this object.
    const FACTORS: Factorization;
}

impl Factorization {
    /// Creates a new factorization from the given prime powers.
    pub const fn new(prime_powers: &'static [(u128, usize)]) -> Factorization {
        let mut value = 1;
        let mut i = 0;
        while i < prime_powers.len() {
            value *= intpow::<0>(prime_powers[i].0, prime_powers[i].1 as u128);
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
    pub fn maximal_divisors<const L: usize>(&self, l: u128) -> impl Iterator<Item = [usize; L]> + '_ {
        FactorStream::new(&self.prime_powers, l, true)
    }

    /// True if there the factorization represents 1.
    /// False otherwise.
    pub const fn is_empty(&self) -> bool {
        self.prime_powers.is_empty()
    }

    /// Returns the prime power factor represented by prime number `i`, $p_i^{t_i}$.
    /// This method will `panic` if `i` is out of bounds.
    pub const fn factor(&self, i: usize) -> u128 {
        intpow::<0>(self.prime_powers[i].0, self.prime_powers[i].1 as u128)
    }

    /// Gets the prime powers as an array.
    /// The first element of each entry is the prime, and the second is the power.
    pub const fn prime_powers(&self) -> &'static [(u128, usize)] {
        self.prime_powers
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub const fn value(&self) -> u128 {
        self.value
    }

    /// Converts an array of powers on the prime factors into an integer.
    pub fn from_powers(&self, ds: &[usize]) -> u128 {
        let mut total = 1;
        for (d, (p, t)) in ds.iter().zip(self.prime_powers) {
            if d > t {
                return 0;
            } else {
                total *= intpow::<0>(*p, *d as u128);
            }
        }
        total
    }

    /// Returns the number of prime factors in the factorization.
    pub const fn len(&self) -> usize {
        self.prime_powers.len()
    }
}

impl Index<usize> for Factorization {
    type Output = (u128, usize);

    fn index(&self, index: usize) -> &(u128, usize) {
        &self.prime_powers[index]
    }
}
