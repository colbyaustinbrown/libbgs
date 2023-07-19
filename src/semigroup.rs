use std::rc::Rc;
use std::fmt::Debug;

pub trait Semigroup {
    type Elem: SemigroupElem;
    fn order(&self) -> u128;
}

pub trait SemigroupElem: Debug + Clone {
    type Group: Semigroup;
    fn one(param: &Rc<Self::Group>) -> Self;
    fn is_one(&self) -> bool;
    fn param(&self) -> &Rc<Self::Group>;
    fn multiply(&mut self, other: &Self);
    fn square(&mut self);

    fn pow(&mut self, mut n: u128) where Self: Sized {
        let mut y = Self::one(self.param());
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

