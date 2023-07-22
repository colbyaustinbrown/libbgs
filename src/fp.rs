use std::rc::Rc;

use crate::semigroup::*;
use crate::util::long_multiply;
use crate::sylow::*;
use crate::util::*;
use crate::factorization::*;

pub type FpStar = Factorization;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FpNumber {
    value: u128,
    group: Rc<FpStar>
}

impl FpStar {
    fn p(&self) -> u128 {
        self.value() + 1
    }

    pub fn from_int(self: &Rc<Self>, value: u128) -> FpNumber {
        FpNumber {
            value, 
            group: Rc::clone(self)
        }
    }
}

impl FpNumber {
    pub fn value(&self) -> u128 {
        self.value
    }
}

impl PartialEq<u128> for FpNumber {
    fn eq(&self, other: &u128) -> bool {
        self.value == *other
    }
}

impl Semigroup for FpStar {
    type Elem = FpNumber;

    fn size(&self) -> &Factorization {
        self
    }

    fn one(self: &Rc<Self>) -> FpNumber {
        self.from_int(1)
    }
}

impl SylowDecomposable for FpStar {
    fn find_sylow_generator(self: &Rc<Self>, i: usize) -> FpNumber {
        match self.size().prime_powers()[i] {
            (2,1) => self.from_int(self.size().value()),
            _ => (1..self.size().value())
                .map(|i| self.from_int(standard_affine_shift(self.p(), i)))
                .find_map(|c| self.is_sylow_generator(&c, i))
                .unwrap()
        }
    }
}

impl SemigroupElem for FpNumber {
    type Group = FpStar;
    fn is_one(&self) -> bool {
        self.value == 1
    }

    fn group(&self) -> &Rc<FpStar> {
        &self.group
    }

    fn multiply(&mut self, other: &FpNumber) {
        self.value = long_multiply(self.value, other.value, self.group.p());
    }

    fn square(&mut self) {
        self.value = long_multiply(self.value, self.value, self.group.p());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factorization::Factorization;
    use crate::sylow::tests::*;

    #[test]
    fn one_is_one() {
        let p = Rc::new(Factorization::new(
            vec![(2,1), (3,1)]
        ));
        let one = p.one();
        assert!(one.is_one());
    }

    #[test]
    fn multiplies() {
        let p = Rc::new(Factorization::new(
            vec![(2,1), (3,1)]
        ));
        let mut x = p.from_int(3);
        let five = p.from_int(5);
        x.multiply(&five);
        assert_eq!(1, x.value);
        assert_eq!(5, five.value);
    }

    #[test]
    fn squares() {
        let p = Rc::new(Factorization::new(
            vec![(2,1), (3,1)]
        ));
        let mut x = p.from_int(3);
        x.square();
        assert_eq!(2, x.value);

        let mut x = p.from_int(2);
        x.square();
        assert_eq!(4, x.value);
    }

    #[test]
    fn powers_up() {
        let p = Rc::new(Factorization::new(
            vec![(2,1), (3,1)]
        ));
        let mut x = p.from_int(2);
        x.pow(5);
        assert_eq!(4, x.value());

        let mut x = p.from_int(3);
        x.pow(3);
        assert_eq!(6, x.value());

        let mut x = p.from_int(5);
        x.pow(p.size().value());
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        ));
        let mut x = fp.from_int(3);
        x.pow(fp.size().value());
        assert!(x.is_one());
    }

    #[test]
    fn sylow_one_is_one() {
        let fp = Rc::new(Factorization::new(
            vec![(2,2), (3,1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        let one = g.one();
        assert!(one.is_one());
    }

    #[test]
    fn sylow_finds_generators() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 2), (5, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().factor(i);
            test_is_generator_small::<FpStar>(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().prime_powers()[i];
            test_is_generator_big::<FpStar>(gen, d);
        }
    }

    #[test]
    fn sylow_order() {
        let fp = Rc::new(Factorization::new(
            vec![(2,2), (3,1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 1..13 {
            let mut x = SylowElem {
                group: Rc::clone(&g),
                coords: vec![i % 4, i % 3]
            };
            x.pow(x.order().value());
            assert!(x.is_one());
        }
    }

    #[test]
    fn sylow_order_big() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        let n = 123456789;
        let mut x = SylowElem {
            group: Rc::clone(&g),
            coords: g.size().prime_powers().iter()
                .map(|(p,d)| n % intpow(*p, *d, 0))
                .collect()
        };
        let or = x.order().value();
        x.pow(or);
        assert!(x.is_one());
    }
}

