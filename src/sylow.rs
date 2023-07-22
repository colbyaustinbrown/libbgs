use std::rc::Rc;
use std::fmt;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: Semigroup {
    fn is_sylow_generator(self: &Rc<Self>, candidate: &Self::Elem, d: &(u128, u128)) -> Option<Self::Elem> {
        let pow = self.size().value / intpow(d.0, d.1, 0);
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
pub struct SylowDecomp<C: SylowDecomposable> {
    pub parent: Rc<C>,
    pub generators: Vec<C::Elem>
}

#[derive(Eq,PartialEq)]
pub struct SylowElem<C: SylowDecomposable> {
    pub group: Rc<SylowDecomp<C>>,
    pub coords: Vec<u128>
}

impl<C: SylowDecomposable> SylowDecomp<C> {
    pub fn new(parent: &Rc<C>) -> SylowDecomp<C> {
        let length = parent.size().len();
        let mut gen = vec![parent.one(); length];
        for i in 0..length {
            gen[i] = parent.find_sylow_generator(&parent.size().prime_powers[i]);
        }
        SylowDecomp {
            parent: Rc::clone(parent),
            generators: gen
        }
    }
}

impl<C: SylowDecomposable> Semigroup for SylowDecomp<C> {
    type Elem = SylowElem<C>;
    fn size(&self) -> &Factorization {
        self.parent.size()
    }
    fn one(self: &Rc<Self>) -> SylowElem<C> {
        SylowElem {
            group: Rc::clone(self),
            coords: vec![0; self.size().len()]
        }
    }
}

impl<C: SylowDecomposable> Clone for SylowElem<C> {
    fn clone(&self) -> Self {
        SylowElem {
            group: Rc::clone(&self.group),
            coords: self.coords.clone()
        }
    }
}

impl<C: SylowDecomposable> fmt::Debug for SylowElem<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

impl<C: SylowDecomposable> SemigroupElem for SylowElem<C> {
    type Group = SylowDecomp<C>;

    fn is_one(&self) -> bool {
        self.coords.iter().all(|x| { *x == 0 })
    }

    fn group(&self) -> &Rc<SylowDecomp<C>> {
        &self.group
    }

    fn multiply(&mut self, other: &SylowElem<C>) {
        let len = self.group().generators.len();
        for i in 0..usize::max(len, len) {
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group().size().factor(i);
        }
    }

    fn square(&mut self) {
        for i in 0..self.group().generators.len() {
            self.coords[i] = self.coords[i] * 2 % self.group.size().factor(i);
        }
    }
}

impl<C: SylowDecomposable> SylowElem<C> {
    pub fn to_product(&self) -> C::Elem {
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
        for i in 0..self.group.size().prime_powers.len() {
            let mut x = self.clone();
            for j in 0..self.group.size().prime_powers.len() {
                if j == i { continue; }
                x.pow(self.group.size().factor(j));
            }
            let mut r = 0;
            while !x.is_one() {
                x.pow(self.group.size().prime_powers[i].0);
                r += 1;
            }
            prime_powers.push((self.group.size().prime_powers[i].0, r));
        }
        let factors: Vec<u128> = prime_powers.iter()
            .map(|(p, r)| intpow(*p, *r, 0))
            .collect();
        let value: u128 = factors.iter().product();
        Factorization {
            value,
            prime_powers
        }
    }
}

impl<C: SylowDecomposable> SylowDecomposable for SylowDecomp<C> {
    fn find_sylow_generator(self: &Rc<Self>, d: &(u128, u128)) -> Self::Elem {
        let mut coords = vec![0 ; self.generators.len()];
        for i in 0..self.generators.len() {
            if *d == self.size().prime_powers[i] {
                coords[i] = 1;
            }
        }
        SylowElem {
            group: Rc::clone(self),
            coords
        }
    }
}

pub mod tests {
    use super::*;

    // utility method for external tests
    pub fn test_is_generator_small<C: SylowDecomposable> (g: &C::Elem, d: u128) -> bool {
        let mut x = g.clone();
        for _ in 1..d {
            if x.is_one() {return false;}
            x.multiply(g);
        }
        x.is_one()
    }

    // utility method for external tests
    pub fn test_is_generator_big<C: SylowDecomposable>(g: &C::Elem, d: (u128, u128)) {
        let mut x = g.clone();
        for _ in 0..d.1 {
            assert!(!x.is_one());
            x.pow(d.0);
        }
        assert!(x.is_one());
    }
}

