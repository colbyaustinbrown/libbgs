use std::rc::Rc;
use std::fmt;

pub trait Semigroup: Eq {
    type Elem: SemigroupElem<Group = Self>;
    fn size(&self) -> u128;
    fn one(self: &Rc<Self>) -> Self::Elem;
}

pub trait SemigroupElem: Clone + PartialEq + Eq + fmt::Debug {
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

