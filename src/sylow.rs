use std::rc::Rc;

use crate::util::*;
use crate::semigroup::*;
use crate::factorization::*;

pub trait SylowDecomposable: SemigroupElem + Copy {
    type Iter: IntoIterator<Item = Self>;
    fn candidates() -> Self::Iter;
    fn sylow_generator(param: &Rc<Self::Group>, d: (u128, u128)) -> Self 
    where Self: Sized {
        let pred = |candidate: &Self| {
            let power = intpow(d.0, d.1 - 1, 0);
            let mut can = *candidate;
            can.pow(param.order() / (power * d.1));
            !can.is_one()
        };
        Self::candidates().into_iter().find(pred).unwrap()
    }
}

pub struct SylowDecomp<T: SylowDecomposable> {
    length: usize,
    order: Factorization,
    generators: Vec<T>
}

impl<T: SylowDecomposable + Clone> SylowDecomp<T> {
    fn new(param: &Rc<T::Group>, order: Factorization) -> SylowDecomp<T> {
        let length = order.len();
        let mut gen = vec![T::one(&param); length];
        for i in 0..length {
            gen[i] = T::sylow_generator(param, order.primepowers[i]);
        }
        SylowDecomp {
            length: length,
            order,
            generators: gen
        }
    }
}

impl<T: SylowDecomposable> Semigroup for SylowDecomp<T> {
    fn order(&self) -> u128 {
        self.order.value
    }
}

pub struct SylowElem<T: SylowDecomposable> {
    group: Rc<SylowDecomp<T>>,
    coords: Vec<u128>
}

impl<T: SylowDecomposable> SemigroupElem for SylowElem<T> {
    type Group = SylowDecomp<T>;
    fn one(param: &Rc<SylowDecomp<T>>) -> SylowElem<T> {
        SylowElem {
            group: Rc::clone(param),
            coords: vec![0; param.length]
        }
    }

    fn is_one(&self) -> bool {
        self.coords.iter().all(|x| { *x == 0 })
    }

    fn param(&self) -> &Rc<SylowDecomp<T>> {
        &self.group
    }

    fn multiply(&mut self, other: &SylowElem<T>) {
        for i in 0..usize::max(self.group.length, other.group.length) {
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group.order.factors[i];
        }
    }

    fn square(&mut self) {
        for i in 0..self.group.length {
            self.coords[i] = self.coords[i] * 2 % self.group.order.factors[i];
        }
    }
}

impl<T: SylowDecomposable + Copy> SylowElem<T> {
    fn to_product(&self, param: &Rc<T::Group>) -> T {
        let mut x = T::one(param);
        for i in 0..self.group.length {
            if self.coords[i] > 0 {
                let mut y = self.group.generators[i];
                y.pow(self.coords[i]);
                x.multiply(&y);
            }
        }
        x
    }
}

