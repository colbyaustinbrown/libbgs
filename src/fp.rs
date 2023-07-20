use std::rc::Rc;

use crate::semigroup::*;
use crate::util::long_multiply;
use crate::sylow::*;
use crate::util::*;

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

impl SylowDecomposable for Fp {
    fn find_sylow_generator(self: &Rc<Self>, d: &(u128, u128)) -> FpNumber {
        match d {
            (2,1) => (**self - 1, Rc::clone(self)),
            _ => {
                (1..**self)
                    .map(|i| (standard_affine_shift(**self, i), Rc::clone(self)))
                    .find_map(|c| self.is_sylow_generator(&c, d))
                    .unwrap()
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::factorization::Factorization;
    use crate::sylow::tests::*;

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

    #[test]
    fn sylow_one_is_one() {
        let fp = Rc::new(13);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 12,
            factors: vec![4, 3],
            prime_powers: vec![(2,2), (3,1)]
        }));
        let one = g.one();
        assert!(one.is_one());
    }

    #[test]
    fn sylow_finds_generators() {
        let fp = Rc::new(13);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 12,
            factors: vec![4, 3],
            prime_powers: vec![(2,2), (3,1)]
        }));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.order_factors.factors[i];
            test_is_generator_small::<Fp>(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp = Rc::new(1_000_000_000_000_000_124_399);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 1_000_000_000_000_000_124_398,
            factors: vec![2, 7, 13, 841, 43, 705737, 215288719],
            prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        }));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.order_factors.prime_powers[i];
            test_is_generator_big::<Fp>(gen, d);
        }
    }

    #[test]
    fn decomposes() {
        let fp = Rc::new(13);
        let g = Rc::new(SylowDecomp::new(&fp, Factorization {
            value: 12,
            factors: vec![4, 3],
            prime_powers: vec![(2,2), (3,1)]
        }));
        println!("generators are {:?}", g.generators);
        let x = (2, Rc::clone(&fp));
        let decomp = g.decompose(&x);
        println!("{:?}", decomp);
        assert_eq!(x, decomp.to_product());
    }
}

