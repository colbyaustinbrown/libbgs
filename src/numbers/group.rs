use std::fmt;

use crate::numbers::Factorization;

pub trait Group: Eq {
    type Elem: GroupElem<Group = Self>;

    fn size() -> u128;
    fn one() -> Self::Elem;
}

pub trait GroupElem: Clone + PartialEq + Eq + fmt::Debug {
    type Group: Group<Elem = Self>;

    fn is_one(&self) -> bool;
    fn multiply(&self, other: &Self) -> Self;
    fn square(&self) -> Self;

    fn pow(&self, mut n: u128) -> Self {
        let mut y = Self::Group::one();
        let mut res = self.clone();
        while n > 1 {
            // println!("{n} {:?} {:?}", &self, y);
            if n % 2 == 1 {
                y = y.multiply(&res);
            }
            res = res.square();
            n >>= 1;
        }
        res.multiply(&y)
    }

    fn inverse(&self) -> Self {
        let res = self.clone();
        res.pow(Self::Group::size() - 1)
    }

    fn order<const L: usize>(&self, parent_size: &Factorization<L>) -> Factorization<L> {
        let mut prime_powers = [(0,0); L];
        for i in 0..L {
            let mut x = self.clone();
            for j in 0..L {
                if j == i {
                    continue;
                }
                x = x.pow(parent_size.factor(j));
            }

            let mut r = 0;
            while !x.is_one() {
                x = x.pow(parent_size[i].0);
                r += 1;
            }
            prime_powers[i] = (parent_size[i].0, r)
        }
        Factorization::new(prime_powers)
    }
}
