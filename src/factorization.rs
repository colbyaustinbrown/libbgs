use std::ops::Index;

pub use crate::util::*;
pub use crate::semigroup::*;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    value: u128,
    prime_powers: Vec<(u128, u128)>,
}

pub trait Factored: {
    fn factors(&self) -> &Factorization;
}

pub trait FactoredElem<G: Factored>: SemigroupElem 
where Self: SemigroupElem<Group = G> {
   fn order(&self) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..self.group().factors().len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..self.group().factors().len() {
                    if j == i { continue; }
                    x.pow(self.group().factors().factor(j));
                }

                let mut r = 0;
                while !x.is_one() {
                    x.pow(self.group().factors()[i].0);
                    r += 1;
                }
                (self.group().factors()[i].0, r)
            })
            .collect();
        Factorization::new(prime_powers)
        
    }
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

    pub fn len(&self) -> usize {
        self.prime_powers.len()
    }

    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    pub fn value(&self) -> u128 {
        self.value
    }

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

