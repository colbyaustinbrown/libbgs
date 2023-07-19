use std::rc::Rc;
use std::fmt;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: Semigroup {
    fn is_sylow_generator(self: &Rc<Self>, candidate: &Self::Elem, d: &(u128, u128)) -> Option<Self::Elem> {
        let respow = self.order() / intpow(d.0, d.1, 0);
        let checkpow = respow * intpow(d.0, d.1 - 1, 0);
        let mut mut_can = candidate.clone();
        mut_can.pow(checkpow);
        if !mut_can.is_one() { 
            mut_can.pow(respow);
            Some(mut_can) 
        } else { None }
    }
    fn find_sylow_generator(self: &Rc<Self>, d: &(u128, u128)) -> Self::Elem;
}

#[derive(Debug)]
pub struct SylowDecomp<G: SylowDecomposable> {
    pub parent: Rc<G>,
    pub order_factors: Factorization,
    pub generators: Vec<G::Elem>
}

impl<G: SylowDecomposable> SylowDecomp<G> {
    pub fn new(parent: &Rc<G>, order_factors: Factorization) -> SylowDecomp<G> {
        println!("Beginning search for generators.");
        println!("Group of order {:?}.", parent.order());
        let length = order_factors.len();
        let mut gen = vec![parent.one(); length];
        for i in 0..length {
            println!("Searching for generator of order {:?}", order_factors.prime_powers[i]);
            gen[i] = parent.find_sylow_generator(&order_factors.prime_powers[i]);
            println!("Found it! {:?}", gen[i]);
        }
        SylowDecomp {
            parent: Rc::clone(parent),
            order_factors,
            generators: gen
        }
    }
}

impl<G: SylowDecomposable> Semigroup for SylowDecomp<G> {
    type Elem = SylowElem<G>;
    fn order(&self) -> u128 {
        self.parent.order() 
    }
    fn one(self: &Rc<Self>) -> SylowElem<G> {
        SylowElem {
            group: Rc::clone(self),
            coords: vec![0; self.order_factors.len()]
        }
    }
}

pub struct SylowElem<G: SylowDecomposable> {
    group: Rc<SylowDecomp<G>>,
    coords: Vec<u128>
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
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group().order_factors.factors[i];
        }
    }

    fn square(&mut self) {
        for i in 0..self.group().generators.len() {
            self.coords[i] = self.coords[i] * 2 % self.group.order_factors.factors[i];
        }
    }
}

impl<G: SylowDecomposable> SylowElem<G> {
    fn to_product(&self) -> G::Elem {
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
}

