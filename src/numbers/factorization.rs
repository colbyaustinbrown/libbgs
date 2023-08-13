use std::ops::Index;

pub use crate::util::*;
pub use crate::numbers::semigroup::*;

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
   fn order(&self, g: &G) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..g.factors().len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..g.factors().len() {
                    if j == i { continue; }
                    x.pow(g.factors().factor(j), g);
                }

                let mut r = 0;
                while !x.is_one() {
                    x.pow(g.factors()[i].0, g);
                    r += 1;
                }
                (g.factors()[i].0, r)
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

    pub fn maximal_divisors(&self, l: u128) -> impl Iterator<Item = Vec<u128>> + '_ {
        FactorStream::new(&self.prime_powers, l)
    }

    pub fn len(&self) -> usize {
        self.prime_powers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.prime_powers.is_empty()
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

struct FactorStream<'a> {
    source: &'a [(u128, u128)],
    stack: Vec<(usize, Vec<u128>)>,
    limit: u128,
    maximal_only: bool
}

impl<'a> FactorStream<'a> {
    fn new(source: &'a [(u128, u128)], limit: u128) -> FactorStream {
        FactorStream {
            source,
            limit,
            stack: vec![(0, vec![0 ; source.len()])],
            maximal_only: true
        }
    }
}

impl<'a> Iterator for FactorStream<'a> {
    type Item = Vec<u128>;

    fn next(&mut self) -> Option<Vec<u128>> {
        let Some((i, state)) = self.stack.pop() else { return None; };
        // println!("{state:?}");
        let prod: u128 = state.iter()
            .zip(self.source)
            .map(|(d, (p, _))| {
                intpow(*p, *d, 0)
            })
            .product();
        let mut maximal = true;
        for j in i..self.source.len() {
            if state[j] == self.source[j].1 { continue; }
            // Assumption: the primes in a factorization are in increasing order.
            // If not, the break below should be a continue.
            if prod * self.source[j].0 > self.limit { break; } 
            let mut next = state.clone();
            next[j] += 1;
            self.stack.push((j, next));
            maximal = false;
        }
        maximal &= prod * self.source[0].0 > self.limit;
        if self.maximal_only && !maximal { self.next() } 
        else { Some(state) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stream() {
        let facts = Factorization::new(vec![(2,3), (3,2), (5,1)]);
        let stream = FactorStream::new(&facts.prime_powers, 25);
        let mut count = 0;
        for x in stream {
            println!("{x:?}");
            // assert!(x < 25);
            count += 1;
        }
        assert_eq!(count, 4);
    }
}

