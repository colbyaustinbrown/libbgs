use std::rc::Rc;

pub trait Semigroup {
    fn order(&self) -> u128;
}

pub trait SemigroupElem {
    type Group: Semigroup;
    fn one(param: &Rc<Self::Group>) -> Self;
    fn is_one(&self) -> bool;
    fn param(&self) -> &Rc<Self::Group>;
    fn multiply(&mut self, other: &Self);
    fn square(&mut self);

    fn pow(&mut self, mut n: u128) where Self: Sized {
        let mut y = Self::one(self.param());
        while n > 0 {
            if n % 2 == 1 {
                y.multiply(self);
            }
            self.square();
            n >>= 1;
        }
        self.multiply(&y);
    }
}

