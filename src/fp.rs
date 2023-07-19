use std::rc::Rc;

use crate::semigroup::*;
use crate::util::long_multiply;

pub type Fp = u128;
pub type FpNumber = (u128, Rc<Fp>);

impl Semigroup for Fp {
    type Elem = FpNumber;
    fn order(&self) -> u128 {
        *self
    }
    fn one(self: &Rc<Self>) -> FpNumber {
        (1, Rc::clone(self))
    }
}

impl SemigroupElem for FpNumber {
    type Group = Fp;
    fn is_one(&self) -> bool {
        self.0 == 1
    }

    fn group(&self) -> &Rc<Fp> {
        &self.1
    }

    fn multiply(&mut self, other: &FpNumber) {
        self.0 = long_multiply(self.0, other.0, *self.1);
    }

    fn square(&mut self) {
        self.0 = long_multiply(self.0, self.0, *self.1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_is_one() {
        let p = Rc::new(7);
        let one = p.one();
        assert!(one.is_one());
    }

    #[test]
    fn multiplies() {
        let p = Rc::new(7);
        let mut x = (3, Rc::clone(&p));
        let five = (5, Rc::clone(&p));
        x.multiply(&five);
        assert_eq!(1, x.0);
        assert_eq!(5, five.0);
    }

    #[test]
    fn squares() {
        let p = Rc::new(7);
        let mut x = (3, Rc::clone(&p));
        x.square();
        assert_eq!(2, x.0);

        let mut x = (2, Rc::clone(&p));
        x.square();
        assert_eq!(4, x.0);
    }

    #[test]
    fn powers_up() {
        let p = Rc::new(7);
        let mut x = (2, Rc::clone(&p));
        x.pow(5);
        assert_eq!(4, x.0);

        let mut x = (3, Rc::clone(&p));
        x.pow(3);
        assert_eq!(6, x.0);

        let mut x = (5, Rc::clone(&p));
        x.pow(*p - 1);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let p = Rc::new((1 << 61) - 1); // Mersenne prime
        let mut x = (3, Rc::clone(&p));
        x.pow(*p - 1);
        assert!(x.is_one());

        let p = Rc::new(1_000_000_000_000_000_124_399);
        let mut x = (3, Rc::clone(&p));
        x.pow(*p - 1);
        assert!(x.is_one());
    }
}

