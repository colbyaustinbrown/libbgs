use std::rc::Rc;

use crate::semigroup::*;
use crate::util::long_multiply;
use crate::sylow::*;
use crate::util::*;
use crate::factorization::*;
use crate::quad_field::*;
use crate::group::*;

pub type FpStar = Factorization;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FpNum {
    value: u128
}

impl FpStar {
    pub fn p(&self) -> u128 {
        self.value() + 1
    }

    pub fn from_int(&self, value: u128) -> FpNum {
        FpNum {
            value
        }
    }
}

impl Factored for FpStar {
    fn factors(&self) -> &Factorization { self }
}

impl FpNum {
    pub fn value(&self) -> u128 {
        self.value
    }

    pub fn int_sqrt(&self, fp: &FpStar) -> Option<FpNum> {
        // Cipolla's algorithm

        let p = fp.p();
        let l = legendre(self.value, fp.p());
        if l == 0 { 
            return Some(fp.from_int(0));
        } else if l == p - 1{
            return None;
        }
        let mut i = 1;
        let (a,r) = loop {
            let a = standard_affine_shift(p, i);
            let r = (intpow(a, 2, p) + p - self.value) % p;
            if legendre(r, p) == (p - 1) { break (a,r); }
            i += 1;
        };
        let fp2 = Rc::new(QuadFieldBare::new(p, r));
        let mut x = fp2.from_ints(a, 1);
        x.pow((p + 1) / 2, &fp2);
        Some(fp.from_int(x.a0))
    }
}

impl PartialEq<u128> for FpNum {
    fn eq(&self, other: &u128) -> bool {
        self.value == *other
    }
}

impl Semigroup for FpStar {
    type Elem = FpNum;

    fn size(&self) -> u128 {
        self.factors().value()
    }

    fn one(&self) -> FpNum {
        self.from_int(1)
    }
}

impl Group for FpStar {}

impl SylowDecomposable for FpStar {
    fn find_sylow_generator(&self, i: usize) -> FpNum {
        match self.factors()[i] {
            (2,1) => self.from_int(self.size()),
            _ => (1..self.size())
                .map(|i| self.from_int(standard_affine_shift(self.p(), i)))
                .find_map(|c| self.is_sylow_generator(&c, i))
                .unwrap()
        }
    }
}

impl SemigroupElem for FpNum {
    type Group = FpStar;
    fn is_one(&self, _: &FpStar) -> bool {
        self.value == 1
    }

    fn multiply(&mut self, other: &FpNum, g: &FpStar) {
        self.value = long_multiply(self.value, other.value, g.p());
    }

    fn square(&mut self, g: &FpStar) {
        self.value = long_multiply(self.value, self.value, g.p());
    }
}

impl FactoredElem<FpStar> for FpNum {}

impl GroupElem for FpNum {
    fn invert(&mut self, g: &FpStar) {
        self.pow(g.p() - 2, g); 
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factorization::Factorization;
    use crate::sylow::tests::*;

    #[test]
    fn one_is_one() {
        let p = Factorization::new(vec![(2,1), (3,1)]);
        let one = p.one();
        assert!(one.is_one(&p));
    }

    #[test]
    fn multiplies() {
        let p = Factorization::new(vec![(2,1), (3,1)]);
        let mut x = p.from_int(3);
        assert_eq!(3, x.value);
        let five = p.from_int(5);
        x.multiply(&five, &p);
        assert_eq!(1, x.value);
        assert_eq!(5, five.value);
    }

    #[test]
    fn squares() {
        let p = Factorization::new(vec![(2,1), (3,1)]);
        let mut x = p.from_int(3);
        x.square(&p);
        assert_eq!(2, x.value);

        let mut x = p.from_int(2);
        x.square(&p);
        assert_eq!(4, x.value);
    }

    #[test]
    fn powers_up() {
        let p = Factorization::new(vec![(2,1), (3,1)]);
        let mut x = p.from_int(2);
        x.pow(5, &p);
        assert_eq!(4, x.value());

        let mut x = p.from_int(3);
        x.pow(3, &p);
        assert_eq!(6, x.value());

        let mut x = p.from_int(5);
        x.pow(p.factors().value(), &p);
        assert!(x.is_one(&p));
    }

    #[test]
    fn powers_up_big() {
        let p = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let mut x = p.from_int(3);
        x.pow(p.value(), &p);
        assert!(x.is_one(&p));
    }

    #[test]
    fn sylow_one_is_one() {
        let p = Factorization::new(vec![(2,2), (3,1)]);
        let g = SylowDecomp::new(&p);
        let one = g.one();
        assert!(one.is_one(&g));
    }

    #[test]
    fn sylow_finds_generators() {
        let p = Factorization::new(vec![(2, 2), (5, 1)]);
        let g = SylowDecomp::new(&p);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.parent.factor(i);
            test_is_generator_small::<FpStar>(gen, d, &p);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let p = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&p);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big::<FpStar>(gen, d, &p);
        }
    }

    #[test]
    fn sylow_order() {
        let fp = Factorization::new(vec![(2,2), (3,1)]);
        let g = SylowDecomp::new(&fp);
        for i in 1..13 {
            let mut x = SylowElem::new(vec![i % 4, i % 3]);
            x.pow(x.order(&g).value(), &g);
            assert!(x.is_one(&g));
        }
    }

    #[test]
    fn sylow_order_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&fp);
        let n = 123456789;
        let mut x = SylowElem::new(
            g.factors().as_array().iter()
                .map(|(p,d)| n % intpow(*p, *d, 0))
                .collect()
        );
        let or = x.order(&g).value();
        x.pow(or, &g);
        assert!(x.is_one(&g));
    }

    #[test]
    fn calculates_square_roots() {
        let p = Factorization::new(vec![(2, 2), (3, 1)]);
        let mut nonresidues = 0;
        for x in (1..13).map(|i| p.from_int(i)) {
            match x.int_sqrt(&p) {
                None => { nonresidues += 1; }
                Some(mut y) => {
                    y.pow(2, &p);
                    assert_eq!(x, y);
                }
            }
        }
        assert_eq!(nonresidues, 6);
    }

    #[test]
    fn inverts() {
        let p = Factorization::new(vec![(2, 2), (3, 1)]);
        for i in 2..13 {
            let mut x = p.from_int(i);
            let y = x.clone();
            x.invert(&p);
            assert!(!x.is_one(&p));
            x.multiply(&y, &p);
            assert!(x.is_one(&p));
        }
    }
}

