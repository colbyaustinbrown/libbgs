use crate::util::*;
use crate::semigroup::*;

pub trait Cyclic {
    fn of_order(n: i128) -> Self;
}

pub struct SylowDecomp<T: Cyclic> {
    length: usize,
    orders: Vec<u128>,
    generators: Vec<Option<T>>
}

pub struct SylowElem<'a, T: Cyclic> {
    group: &'a SylowDecomp<T>,
    coords: Vec<u128>
}

impl<'a,T: Cyclic> ParameterizedSemigroup<&'a SylowDecomp<T>> for SylowElem<'a, T> {
    fn one(param: &'a SylowDecomp<T>) -> SylowElem<'a, T> {
        SylowElem {
            group: param,
            coords: vec![0; param.length]
        }
    }

    fn param(&self) -> &'a SylowDecomp<T> {
        self.group
    }

    fn multiply(&mut self, other: &SylowElem<'a, T>) {
        for i in 0..usize::max(self.group.length, other.group.length) {
            self.coords[i] = (self.coords[i] + other.coords[i]) % self.group.orders[i];
        }
    }

    fn square(&mut self) {
        for i in 0..self.group.length {
            self.coords[i] = self.coords[i] * 2 % self.group.orders[i];
        }
    }
}


