use std::fmt;
use std::marker::PhantomData;

use crate::util::*;
pub use crate::numbers::semigroup::*;
use crate::numbers::factorization::*;
use crate::numbers::group::*;

pub trait SylowDecomposable: Semigroup + Factored {
    fn find_sylow_generator(&self, i: usize) -> Self::Elem;

    fn is_sylow_generator(&self, candidate: &Self::Elem, i: usize) -> Option<Self::Elem> {
        let d = self.factors()[i];
        let pow = self.size() / intpow(d.0, d.1, 0);
        let mut res = candidate.clone();
        res.pow(pow, self);
        let mut check = res.clone();
        check.pow(intpow(d.0, d.1 - 1, 0), self);
        if check.is_one(self) { None } else { Some(res) }
    }
}

#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<'a, C: SylowDecomposable> {
    pub parent: &'a C,
    pub generators: Vec<C::Elem>
}

#[derive(Eq,PartialEq)]
pub struct SylowElem<'a, C: SylowDecomposable> {
    _group: PhantomData<fn(&'a C) -> &'a C>,
    pub coords: Vec<u128>
}

impl<'a, C: SylowDecomposable> SylowDecomp<'a, C> {
    pub fn new(parent: &C) -> SylowDecomp<C> {
        let length = parent.factors().len();
        let generators = (0..length)
            .map(|i| parent.find_sylow_generator(i))
            .collect();
        SylowDecomp {
            parent,
            generators
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }
    
    fn len(&self) -> usize {
        self.generators.len()
    }
}

impl<'a, C: SylowDecomposable> Semigroup for SylowDecomp<'a, C> {
    type Elem = SylowElem<'a, C>;

    fn one(&self) -> SylowElem<'a, C> {
        SylowElem {
            _group: PhantomData,
            coords: vec![0; self.factors().len()]
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }
}

impl<'a, C: SylowDecomposable> Group for SylowDecomp<'a, C> {}

impl<'a, C: SylowDecomposable> Factored for SylowDecomp<'a, C> {
    fn factors(&self) -> &Factorization {
        self.parent.factors()
    }
}

impl<'a, C: SylowDecomposable> SylowDecomposable for SylowDecomp<'a, C> {
    fn find_sylow_generator(&self, i: usize) -> Self::Elem {
        let mut coords = vec![0 ; self.factors().len()];
        coords[i] = 1;
        SylowElem {
            _group: PhantomData,
            coords
        }
    }
}

impl<'a, C: SylowDecomposable> SylowElem<'a, C> {
    pub fn new(coords: Vec<u128>) -> SylowElem<'a, C> {
        SylowElem {
            _group: PhantomData,
            coords
        }
    }
    pub fn to_product(&self, g: &SylowDecomp<C>) -> C::Elem {
        (0..g.len())
            .filter(|i| self.coords[*i] > 0)
            .fold(g.parent.one(), |mut x, i| {
                let mut y = g.generators[i].clone();
                y.pow(self.coords[i], g.parent);
                x.multiply(&y, g.parent);
                x
            })
    }

    pub fn order(&self, g: &'a SylowDecomp<C>) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..g.len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..g.len() {
                    if j == i { continue; }
                    x.pow(g.factors().factor(j), g);
                }

                let mut r = 0;
                while !x.is_one(g) {
                    x.pow(g.factors()[i].0, g);
                    r += 1;
                }
                (g.factors()[i].0, r)
            })
            .collect();
        Factorization::new(prime_powers)
    }
}

impl<'a, C> SemigroupElem for SylowElem<'a, C> 
where C: SylowDecomposable + 'a {
    type Group = SylowDecomp<'a, C>;

    fn is_one(&self, _: &SylowDecomp<C>) -> bool {
        self.coords.iter().all(|x| { *x == 0 })
    }

    fn multiply(&mut self, other: &SylowElem<C>, g: &SylowDecomp<C>) {
        let len = g.generators.len();
        for i in 0..usize::max(len, len) {
            self.coords[i] = (self.coords[i] + other.coords[i]) % g.factors().factor(i);
        }
    }

    fn square(&mut self, g: &SylowDecomp<C>) {
        for i in 0..g.generators.len() {
            self.coords[i] = self.coords[i] * 2 % g.factors().factor(i);
        }
    }
}

impl<'a, C: SylowDecomposable> GroupElem for SylowElem<'a, C> {
    fn invert(&mut self, g: &Self::Group) {
        for i in 0..self.coords.len() {
            self.coords[i] = g.factors().factor(i) - self.coords[i];
        }
    }
}

impl<'a, C: SylowDecomposable> Clone for SylowElem<'a, C> {
    fn clone(&self) -> Self {
        SylowElem {
            _group: PhantomData,
            coords: self.coords.clone()
        }
    }
}

impl<'a, C: SylowDecomposable> fmt::Debug for SylowElem<'a, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

pub mod tests {
    use super::*;

    // utility method for external tests
    pub fn test_is_generator_small<C: SylowDecomposable> (x: &C::Elem, d: u128, g: &C) -> bool {
        let mut y = x.clone();
        for _ in 1..d {
            if y.is_one(g) {return false;}
            y.multiply(x, g);
        }
        y.is_one(g)
    }

    // utility method for external tests
    pub fn test_is_generator_big<C: SylowDecomposable>(x: &C::Elem, d: (u128, u128), g: &C) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(!y.is_one(g));
            y.pow(d.0, g);
        }
        y.pow(d.0, g);
        assert!(y.is_one(g));
    }
}

