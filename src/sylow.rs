use std::rc::Rc;
use std::fmt;

use crate::util::*;
pub use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: Semigroup + Factorized {
    fn find_sylow_generator(self: &Rc<Self>, i: usize) -> Self::Elem;

    fn is_sylow_generator(self: &Rc<Self>, candidate: &Self::Elem, i: usize) -> Option<Self::Elem> {
        let d = self.factors()[i];
        let pow = self.size() / intpow(d.0, d.1, 0);
        let mut res = candidate.clone();
        res.pow(pow);
        let mut check = res.clone();
        check.pow(intpow(d.0, d.1 - 1, 0));
        if check.is_one() { None } else { Some(res) }
    }
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
        let length = parent.factors().len();
        let gen = (0..length)
            .map(|i| {
                let res = parent.find_sylow_generator(i);
                res
            })
            .collect();
        SylowDecomp {
            parent: Rc::clone(parent),
            generators: gen
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }
}

impl<C: SylowDecomposable> Semigroup for SylowDecomp<C> {
    type Elem = SylowElem<C>;

    fn one(self: &Rc<Self>) -> SylowElem<C> {
        SylowElem {
            group: Rc::clone(self),
            coords: vec![0; self.factors().len()]
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }
}

impl<C: SylowDecomposable> Factorized for SylowDecomp<C> {
    fn factors(&self) -> &Factorization {
        self.parent.factors()
    }
}

impl<C: SylowDecomposable> SylowDecomposable for SylowDecomp<C> {
    fn find_sylow_generator(self: &Rc<Self>, i: usize) -> Self::Elem {
        let mut coords = vec![0 ; self.factors().len()];
        coords[i] = 1;
        SylowElem {
            group: Rc::clone(self),
            coords
        }
    }
}

impl<C: SylowDecomposable> SylowElem<C> {
    pub fn to_product(&self) -> C::Elem {
        (0..self.len())
            .filter(|i| self.coords[*i] > 0)
            .fold(self.group().parent.one(), |mut x, i| {
                let mut y = self.group.generators[i].clone();
                y.pow(self.coords[i]);
                x.multiply(&y);
                x
            })
    }

    pub fn order(&self) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..self.len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..self.len() {
                    if j == i { continue; }
                    x.pow(self.group.factors().factor(j));
                }

                let mut r = 0;
                while !x.is_one() {
                    x.pow(self.group.factors()[i].0);
                    r += 1;
                }
                (self.group.factors()[i].0, r)
            })
            .collect();
        Factorization::new(prime_powers)
    }

    pub fn len(&self) -> usize {
        self.group().factors().len()
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
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group().factors().factor(i);
        }
    }

    fn square(&mut self) {
        for i in 0..self.group().generators.len() {
            self.coords[i] = self.coords[i] * 2 % self.group.factors().factor(i);
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

