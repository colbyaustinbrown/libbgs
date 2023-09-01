use std::ops::Index;

use crate::numbers::FactorStream;
use crate::util::intpow;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization<const L: usize> {
    value: u128,
    prime_powers: [(u128, u128); L],
}

pub trait Factored<S, const L: usize> {
    const FACTORS: Factorization<L>;
}

impl<const L: usize> Factorization<L> {
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

    pub fn maximal_divisors(&self, l: u128) -> impl Iterator<Item = [u128; L]> + '_ {
        FactorStream::new(&self.prime_powers, l)
    }

    pub fn is_empty(&self) -> bool {
        self.prime_powers.is_empty()
    }

    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    pub fn prime_powers(&self) -> &[(u128, u128); L] {
        &self.prime_powers
    }

    pub fn value(&self) -> u128 {
        self.value
    }

    pub fn as_array(&self) -> &[(u128, u128)] {
        &self.prime_powers
    }
}

impl<const L: usize> Index<usize> for Factorization<L> {
    type Output = (u128, u128);

    fn index(&self, index: usize) -> &(u128, u128) {
        &self.prime_powers[index]
    }
}
