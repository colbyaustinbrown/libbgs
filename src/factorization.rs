use crate::util::*;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    value: u128,
    prime_powers: Vec<(u128, u128)>,
}

impl Factorization {
    pub fn new(prime_powers: Vec<(u128, u128)>) -> Factorization {
        Factorization {
            value: prime_powers.iter()
                .map(|(p, r)| intpow(*p, *r, 0))
                .product(),
            prime_powers
        }
    }

    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    pub fn value(&self) -> u128 {
        self.value
    }

    pub fn prime_powers(&self) -> &Vec<(u128, u128)> {
        &self.prime_powers
    }

    pub fn len(&self) -> usize {
        self.prime_powers.len()
    }
}

