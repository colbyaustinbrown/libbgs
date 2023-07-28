use std::fmt;

pub trait Semigroup: Eq {
    type Elem: SemigroupElem<Group = Self>;
    fn size(&self) -> u128;
    fn one(&self) -> Self::Elem;
}

pub trait SemigroupElem: Clone + PartialEq + Eq + fmt::Debug {
    type Group: Semigroup<Elem = Self>;
    fn is_one(&self, g: &Self::Group) -> bool;
    fn multiply(&mut self, other: &Self, g: &Self::Group);
    fn square(&mut self, g: &Self::Group);

    fn pow(&mut self, mut n: u128, g: &Self::Group) where Self: Sized {
        let mut y = Self::Group::one(g);
        while n > 1 {
            // println!("{n} {:?} {:?}", &self, y);
            if n % 2 == 1 {
                y.multiply(self, g);
            }
            self.square(g);
            n >>= 1;
        }
        self.multiply(&y, g);
    }
}

