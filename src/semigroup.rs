use std::rc::Rc;
use std::fmt;

use crate::factorization::*;

pub trait Semigroup: Eq {
    type Elem: SemigroupElem<Group = Self>;
    fn size(&self) -> &Factorization;
    fn one(self: &Rc<Self>) -> Self::Elem;

    fn factors(&self) -> &Vec<(u128, u128)> {
        self.size().prime_powers()
    }
}

pub trait SemigroupElem: Clone + Eq + fmt::Debug {
    type Group: Semigroup<Elem = Self>;
    fn is_one(&self) -> bool;
    fn group(&self) -> &Rc<Self::Group>;
    fn multiply(&mut self, other: &Self);
    fn square(&mut self);

    fn pow(&mut self, mut n: u128) where Self: Sized {
        let mut y = Self::Group::one(self.group());
        while n > 1 {
            // println!("{n} {:?} {:?}", &self, y);
            if n % 2 == 1 {
                y.multiply(self);
            }
            self.square();
            n >>= 1;
        }
        // println!("{:?} {:?}", &self, y);
        self.multiply(&y);
        // println!("{:?} {:?}", &self, y);
    }
}

