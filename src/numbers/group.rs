use std::fmt;

use crate::numbers::factorization::*;

pub trait Group: Eq {
    type Elem: GroupElem<Group = Self>;

    fn size(&self) -> u128;
    fn one(&self) -> Self::Elem;
}

pub trait GroupElem: Clone + PartialEq + Eq + fmt::Debug {
    type Group: Group<Elem = Self>;

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

    fn invert(&self, g: &Self::Group) -> Self {
        let res = self.clone();
        res.pow(g.size() - 1, g)
    }

    fn order(&self, parent: &Self::Group, parent_size: &Factorization) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..parent_size.len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..parent_size.len() {
                    if j == i {
                        continue;
                    }
                    x = x.pow(parent_size.factor(j), parent);
                }

                let mut r = 0;
                while !x.is_one(parent) {
                    x = x.pow(parent_size[i].0, parent);
                    r += 1;
                }
                (parent_size[i].0, r)
            })
            .collect();
        Factorization::new(prime_powers)
    }
}
