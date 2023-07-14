use crate::quad_field::*;

pub trait ParameterizedSemigroup<T> {
    fn one(param: T) -> Self;
    fn param(&self) -> T;
    fn multiply(&mut self, other: &Self);
    fn square(&mut self);

    fn pow(&mut self, mut n: i128) where Self: Sized {
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

