use std::ops::Index;

pub use crate::numbers::group::*;
pub use crate::util::*;

use crate::numbers::factor_stream::*;

/// A prime power decomposition of a positive integer.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    value: u128,
    prime_powers: Vec<(u128, u128)>,
}

/// Types that have a size or order which can be expressed as a product of prime powers.
pub trait Factored {
    /// The prime factorization of this object.
    fn factors(&self) -> &Factorization;
}

impl Factorization {
    /// Returns the factorization given by the prime powers.
    /// Each element of the vector is the prime and its power, respectively.
    /// That is, if the argument is `[(p_1, t_1), .., (p_n, t_n)]`, then this `Factorization`
    /// represents the number $\prod_{i = 1}^n p_i^{t_i}$.
    pub fn new(prime_powers: Vec<(u128, u128)>) -> Factorization {
        Factorization {
            value: prime_powers
                .iter()
                .map(|(p, r)| intpow(*p, *r, 0))
                .product(),
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
    pub fn maximal_divisors(&self, l: u128) -> impl Iterator<Item = Vec<u128>> + '_ {
        FactorStream::new(&self.prime_powers, l)
    }

    /// Returns the number of unique primes which divide the number.
    pub fn len(&self) -> usize {
        self.prime_powers.len()
    }

    /// True if this `Factorization` represents 1, that is, there are no prime divisors of the
    /// number.
    pub fn is_empty(&self) -> bool {
        self.prime_powers.is_empty()
    }

    /// Returns the prime power factor represented by prime number `i`, $p_i^{t_i}$.
    /// This method will `panic` if `i` is out of bounds.
    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub fn value(&self) -> u128 {
        self.value
    }

    /// Returns the prime power decomposition as an array.
    pub fn as_array(&self) -> &[(u128, u128)] {
        &self.prime_powers
    }
}

impl Index<usize> for Factorization {
    type Output = (u128, u128);

    fn index(&self, index: usize) -> &(u128, u128) {
        &self.prime_powers[index]
    }
}
