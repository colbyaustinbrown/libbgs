use std::rc::Rc;
use std::fmt;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;
use crate::fp::*;
use crate::quad_field::*;

pub trait SylowDecomposable: Semigroup {
    type Iter: IntoIterator<Item = Self::Elem>;
    fn candidates(self: &Rc<Self>) -> Self::Iter;
    fn find_sylow_generator(self: &Rc<Self>, d: (u128, u128)) -> Self::Elem {
        let respow = self.order() / intpow(d.0, d.1, 0);
        let checkpow = respow * intpow(d.0, d.1 - 1, 0);
        let pred = |candidate: &Self::Elem| {
            let mut can = candidate.clone();
            can.pow(checkpow);
            !can.is_one()
        };
        let mut res = self.candidates()
            .into_iter()
            .find(pred)
            .unwrap();
        res.pow(respow);
        res
    }
}

#[derive(Debug)]
pub struct SylowDecomp<G: SylowDecomposable> {
    parent: Rc<G>,
    order_factors: Factorization,
    generators: Vec<G::Elem>
}

impl<G: SylowDecomposable> SylowDecomp<G> {
    fn new(parent: &Rc<G>, order_factors: Factorization) -> SylowDecomp<G> {
        let length = order_factors.len();
        let mut gen = vec![parent.one(); length];
        for i in 0..length {
            gen[i] = parent.find_sylow_generator(order_factors.prime_powers[i]);
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

impl SylowDecomposable for Fp {
    type Iter = Vec<FpNumber>;
    fn candidates(self: &Rc<Self>) -> Self::Iter {
        let mut res = Vec::new();
        for i in 2..**self {
            res.push((i, Rc::clone(self)));
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_is_one() {
        let fp = Rc::new(13);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 12,
            factors: vec![4, 3],
            prime_powers: vec![(2,2), (3,1)]
        }));
        let one = g.one();
        assert!(one.is_one());
    }

    #[test]
    fn finds_generators() {
        let fp = Rc::new(13);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 12,
            factors: vec![4, 3],
            prime_powers: vec![(2,2), (3,1)]
        }));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let mut x = gen.clone();
            for _ in 1..g.order_factors.factors[i] {
                println!("{:?}", x);
                assert!(!x.is_one());
                x.multiply(gen);
            }
            println!("{:?}", x);
            assert!(x.is_one());
        }
    }

    #[test]
    fn finds_generators_big() {
        let pplusone = Factorization {
            value: 1_000_000_000_000_000_124_400,
            factors: vec![16, 3, 25, 121, 17, 19, 23, 97, 757, 1453, 8689],
            prime_powers: vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
        };
        let fp = Rc::new(QuadFieldExt::new(
            Factorization {
                value: 1_000_000_000_000_000_124_398,
                factors: vec![2, 7, 13, 841, 43, 705737, 215288719],
                prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            },
            pplusone
        ));
        let g = Rc::new(SylowDecomp::new(&fp, pplusone));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let mut x = gen.clone();
            for _ in 1..g.order_factors.factors[i] {
                println!("{:?}", x);
                assert!(!x.is_one());
                x.multiply(gen);
            }
            println!("{:?}", x);
            assert!(x.is_one());
        }
    }
}

