use std::fmt;

pub trait Semigroup: Eq {
    type Elem: SemigroupElem<Group = Self>;

    fn size(&self) -> u128;
    fn one(&self) -> Self::Elem;
}

pub trait SemigroupElem: Clone + PartialEq + Eq + fmt::Debug {
    type Group: Semigroup<Elem = Self>;

    fn is_one(&self, g: &Self::Group) -> bool;
    fn multiply(&self, other: &Self, g: &Self::Group) -> Self;
    fn square(&self, g: &Self::Group) -> Self;

    fn pow(&self, mut n: u128, g: &Self::Group) -> Self {
        let mut y = Self::Group::one(g);
        let mut res = self.clone();
        while n > 1 {
            // println!("{n} {:?} {:?}", &self, y);
            if n % 2 == 1 {
                y = y.multiply(&res, g);
            }
            res = res.square(g);
            n >>= 1;
        }
        res.multiply(&y, g)
    }
}
