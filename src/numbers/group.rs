use std::fmt;

use crate::numbers::Factorization;

/// Types that represent groups, satisfying the group axioms:
/// * Has an associative binary operator (multiplication)
/// * Has an identity
/// * Every element has an inverse
/// Only the second axiom is enforced at the `Group` trait level, and associativity is never
/// enforced.
/// 
/// Mathematically, this trait can be understood as defining the category $\text{Grp}$, and a struct which
/// implements this trait is a subcategory of $\text{Grp}$. An instance of one of these `struct`s is
/// then a group itself.
pub trait Group: Eq {
    /// The type of the members of this group.
    type Elem: GroupElem<Group = Self>;

    /// Returns the number of elements of this group.
    fn size(&self) -> u128;

    /// Returns the multiplicative identity of the group.
    fn one(&self) -> Self::Elem;
}

/// Types that represent the elements of a group.
///
/// Most methods here are "context dependent", in that they require the group to which the elements
/// belong to be passed as a parameter. Since a `struct` implementing `Group` is more accurately
/// called a subcategory of $\text{Grp}$, it is possible to pass a `Group` to which this
/// `GroupElem` does not actually belong. It is your responsibility to ensure this does not happen.
pub trait GroupElem: Clone + PartialEq + Eq + fmt::Debug {
    /// The group to which these elements belong.
    type Group: Group<Elem = Self>;

    /// True if this element is the multiplicative identity; false otherwise.
    fn is_one(&self, g: &Self::Group) -> bool;

    /// Returns the product of two elements under the group binary operator.
    /// You must guarantee that this operation is associative:
    /// `a.multiply(b.multiply(c, &g), &g) == a.multiply(b, &g).multiply(c, &g)`.
    fn multiply(&self, other: &Self, g: &Self::Group) -> Self;

    /// Returns this element multiplied by itself.
    fn square(&self, g: &Self::Group) -> Self {
        self.multiply(self, g)
    }

    /// Raises this element to the power of `n`.
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

    /// Returns the multiplicative inverse of this element.
    fn inverse(&self, g: &Self::Group) -> Self {
        let res = self.clone();
        res.pow(g.size() - 1, g)
    }

    /// Returns the order of this element, that is, the smallest positive power `p` for which
    /// `a.pow(p, &g).is_one(&g)` returns True.
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
