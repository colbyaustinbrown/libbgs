use std::rc::Rc;
use std::fmt;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: Semigroup {
    fn is_sylow_generator(self: &Rc<Self>, candidate: &Self::Elem, d: &(u128, u128)) -> Option<Self::Elem> {
        let pow = self.order() / intpow(d.0, d.1, 0);
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
    pub order_factors: Factorization,
    pub generators: Vec<G::Elem>
}

impl<G: SylowDecomposable> SylowDecomp<G> {
    pub fn new(parent: &Rc<G>, order_factors: Factorization) -> SylowDecomp<G> {
        let length = order_factors.len();
        let mut gen = vec![parent.one(); length];
        for i in 0..length {
            gen[i] = parent.find_sylow_generator(&order_factors.prime_powers[i]);
        }
        SylowDecomp {
            parent: Rc::clone(parent),
            order_factors,
            generators: gen
        }
    }

    pub fn decompose(self: &Rc<SylowDecomp<G>>, x: &G::Elem) -> SylowElem<G> {
        let mut coords = Vec::new();
        'outer: for i in 0..self.generators.len() {
            let (p,d) = self.order_factors.prime_powers[i];
            let gen = &self.generators[i];

            let precomp_count = if d == 1 {1} else {intpow(p, (d + 1) / 2, 0)};
            let mut precomp = Vec::new();
            let mut a = gen.clone();
            let mut skip = gen.clone();
            skip.pow(precomp_count);
            precomp.push(a.clone());
            for _ in 1..precomp_count {
                a.multiply(&skip);
                precomp.push(a.clone());
            }

            println!("precomputation: {:?}", precomp);

            let mut a = x.clone();
            let killpow = self.order_factors.value / intpow(p, d, 0);
            a.pow(killpow);
            println!("a is {:?}", a);
            for i in 0.. {
                if a.is_one() {
                    coords.push(i);
                    println!("coord is {}", i);
                    continue 'outer;
                }
                for j in 0..precomp.len() {
                    if a == precomp[j] { 
                        println!("j is {}", j);
                        coords.push(precomp_count * (j as u128) - i);     
                        continue 'outer;
                    }
                }
                a.multiply(gen);
            }
        }
        SylowElem {
            group: Rc::clone(self),
            coords
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

#[derive(Eq,PartialEq)]
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
        println!("testing if {:?} has order {:?}", &x, d);
        for _ in 0..d.1 {
            assert!(!x.is_one());
            x.pow(d.0);
            println!("now {:?}", &x);
        }
        assert!(x.is_one());
    }
}

