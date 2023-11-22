use std::ops::Index;

use crate::streams::DivisorStream;
use libbgs_util::intpow;

/// When called with phantom type marker `Ph` and a list of integers, each integer `P` is turned
/// into an implementation of `Factor<Ph> for FpNum<P>` and `Factor<Ph> for QuadNum<P>`.
#[macro_export]
macro_rules! impl_factors {
    ($mrk:ident, $($n:literal),+ $(,)?) => {$(
        impl Factor<$mrk> for FpNum<$n> {
            const FACTORS: Factorization = Factorization::new(make_factor!({$n - 1}));
        }
        impl Factor<$mrk> for QuadNum<$n> {
            const FACTORS: Factorization = Factorization::new(make_factor!({$n + 1}));
        }
    )+};
    ($mrk:ident, $start:literal..$end:literal) => {
        primes!(impl_factors, $mrk, $start..$end);
    }
}
pub use impl_factors;

/// A prime power decomposition of a positive integer.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    value: u128,
    factors: &'static [(u128, usize)],
}

/// A placeholder trait for storing the length of a prime factorization; required to be separate
/// from the Factor trait in order for the trait solver to avoid cycles.
pub trait Length<S> {
    /// The length of a prime factorization, i.e., the number of unique prime factors.
    const LENGTH: usize;
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
    /// The prime factorization of this object.
    const FACTORS: Factorization;
}

impl<S, T: Factor<S>> Length<S> for T {
    const LENGTH: usize = T::FACTORS.len();
}

impl Factorization {
    /// Creates a new factorization from the given prime powers.
    pub const fn new(factors: &'static [(u128, usize)]) -> Factorization {
        let mut value = 1;
        let mut i = 0;
        while i < factors.len() {
            value *= intpow::<0>(factors[i].0, factors[i].1 as u128);
            i += 1;
        }
        Factorization {
            value,
            factors,
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
    pub fn maximal_divisors<const L: usize>(&self, l: u128) -> impl Iterator<Item = [usize; L]> {
        DivisorStream::new(self.factors, l, true).map(|v| v.try_into().unwrap())
    }

    /// True if there the factorization represents 1.
    /// False otherwise.
    pub const fn is_empty(&self) -> bool {
        self.factors.is_empty()
    }

    /// Returns the prime power factor represented by prime number `i`, $p_i^{t_i}$.
    /// This method will `panic` if `i` is out of bounds.
    pub const fn factor(&self, i: usize) -> u128 {
        intpow::<0>(self.factors[i].0, self.factors[i].1 as u128)
    }

    /// Gets the prime powers as an array.
    /// The first element of each entry is the prime, and the second is the power.
    pub const fn factors(&self) -> &'static [(u128, usize)] {
        self.factors
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub const fn value(&self) -> u128 {
        self.value
    }

    /// Converts an array of powers on the prime factors into an integer.
    pub fn from_powers(&self, ds: &[usize]) -> u128 {
        let mut total = 1;
        for (d, (p, t)) in ds.iter().zip(self.factors) {
            if d > t {
                return 0;
            } else {
                total *= intpow::<0>(*p, *d as u128);
            }
        }
        total
    }

    /// Returns the exponents on the factorization.
    pub fn exponents(&self) -> Box<[usize]> {
        self.factors.iter().map(|(_, t)| *t).collect::<Box<[usize]>>()
    }

    /// Returns the number of prime factors in the factorization.
    pub const fn len(&self) -> usize {
        self.factors.len()
    }

    /// Returns $\tau$(`&self`), the number of divisors of this integer.
    /// See OEIS sequence [A000005].
    ///
    /// [A000005]: https://oeis.org/A000005
    pub const fn tau(&self) -> u128 {
        let mut res = 1u128;
        let mut i = 0;
        while i < self.factors.len() {
            res *= (self.factors[i].1 + 1) as u128;
            i += 1;
        }
        res
    }

    /// Returns $\phi$(`&self`), the Euler totient function of this integer.
    /// The totient function $\phi(n)$ is the number of integers $\leq n$ and relatively prime to
    /// $n$. See OEIS sequence [A000010].
    ///
    /// [A000010]: https://oeis.org/A000010
    pub const fn phi(&self) -> u128 {
        let mut res = 1;
        let mut i = 0;
        while i < self.factors.len() {
            res *= intpow::<0>(self.factors[i].0, (self.factors[i].1 - 1) as u128);
            res *= self.factors[i].0 - 1;
            i += 1;
        }
        res
    }
}

impl Index<usize> for Factorization {
    type Output = (u128, usize);

    fn index(&self, index: usize) -> &(u128, usize) {
        &self.factors[index]
    }
}
