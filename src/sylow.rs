use std::rc::Rc;
use std::fmt;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: Semigroup {
    fn is_sylow_generator(self: &Rc<Self>, candidate: &Self::Elem, d: &(u128, u128)) -> Option<Self::Elem> {
        let pow = self.size() / intpow(d.0, d.1, 0);
        let mut res = candidate.clone();
        res.pow(pow);
        let mut check = res.clone();
        check.pow(intpow(d.0, d.1 - 1, 0));
        if !check.is_one() { 
            Some(res) 
        } else { None }
    }
    fn find_sylow_generator(self: &Rc<Self>, d: &(u128, u128)) -> Self::Elem;
}

#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<G: SylowDecomposable> {
    pub parent: Rc<G>,
    pub size: Factorization,
    pub generators: Vec<G::Elem>
}

#[derive(Eq,PartialEq)]
pub struct SylowElem<G: SylowDecomposable> {
    pub group: Rc<SylowDecomp<G>>,
    pub coords: Vec<u128>
}

impl<G: SylowDecomposable> SylowDecomp<G> {
    pub fn new(parent: &Rc<G>, size: Factorization) -> SylowDecomp<G> {
        let length = size.len();
        let mut gen = vec![parent.one(); length];
        for i in 0..length {
            gen[i] = parent.find_sylow_generator(&size.prime_powers[i]);
        }
        SylowDecomp {
            parent: Rc::clone(parent),
            size,
            generators: gen
        }
    }
}

impl<G: SylowDecomposable> Semigroup for SylowDecomp<G> {
    type Elem = SylowElem<G>;
    fn size(&self) -> u128 {
        self.parent.size() 
    }
    fn one(self: &Rc<Self>) -> SylowElem<G> {
        SylowElem {
            group: Rc::clone(self),
            coords: vec![0; self.size.len()]
        }
    }
}

impl<G: SylowDecomposable> Clone for SylowElem<G> {
    fn clone(&self) -> Self {
        SylowElem {
            group: Rc::clone(&self.group),
            coords: self.coords.clone()
        }
    }
}

impl<G: SylowDecomposable> fmt::Debug for SylowElem<G> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

impl<G: SylowDecomposable> SemigroupElem for SylowElem<G> {
    type Group = SylowDecomp<G>;

    fn is_one(&self) -> bool {
        self.coords.iter().all(|x| { *x == 0 })
    }

    fn group(&self) -> &Rc<SylowDecomp<G>> {
        &self.group
    }

    fn multiply(&mut self, other: &SylowElem<G>) {
        let len = self.group().generators.len();
        for i in 0..usize::max(len, len) {
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group().size.factors[i];
        }
    }

    fn square(&mut self) {
        for i in 0..self.group().generators.len() {
            self.coords[i] = self.coords[i] * 2 % self.group.size.factors[i];
        }
    }
}

impl<G: SylowDecomposable> SylowElem<G> {
    pub fn to_product(&self) -> G::Elem {
        let mut x = self.group().parent.one();
        for i in 0..self.group().generators.len() {
            if self.coords[i] > 0 {
                let mut y = self.group.generators[i].clone();
                y.pow(self.coords[i]);
                x.multiply(&y);
            }
        }
        x
    }

    pub fn order(&self) -> Factorization {
        let mut prime_powers = Vec::new(); 
        for i in 0..self.group.size.prime_powers.len() {
            println!("calculating order component of {:?}", self.group.size.prime_powers[i]);
            let mut x = self.clone();
            for j in 0..self.group.size.prime_powers.len() {
                if j == i { continue; }
                x.pow(self.group.size.factors[j]);
            }
            let mut r = 0;
            while !x.is_one() {
                x.pow(self.group.size.prime_powers[i].0);
                r += 1;
            }
            prime_powers.push((self.group.size.prime_powers[i].0, r));
        }
        let factors: Vec<u128> = prime_powers.iter()
            .map(|(p, r)| intpow(*p, *r, 0))
            .collect();
        let value: u128 = factors.iter().product();
        Factorization {
            value,
            factors,
            prime_powers
        }
    }
}

pub mod tests {
    use super::*;

    pub fn test_is_generator_small<G: SylowDecomposable> (g: &G::Elem, d: u128) -> bool {
        let mut x = g.clone();
        for _ in 1..d {
            if x.is_one() {return false;}
            x.multiply(g);
        }
        x.is_one()
    }

    pub fn test_is_generator_big<G: SylowDecomposable>(g: &G::Elem, d: (u128, u128)) {
        let mut x = g.clone();
        for _ in 0..d.1 {
            assert!(!x.is_one());
            x.pow(d.0);
        }
        assert!(x.is_one());
    }
}

