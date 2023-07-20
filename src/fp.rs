use std::rc::Rc;

use crate::semigroup::*;
use crate::util::long_multiply;
use crate::sylow::*;
use crate::util::*;
use crate::factorization::*;

pub type FpStar = Factorization;
pub type FpNumber = (u128, Rc<FpStar>);

impl FpStar {
    fn p(&self) -> u128 {
        self.value + 1
    }
}

impl Semigroup for FpStar {
    type Elem = FpNumber;
    fn size(&self) -> &Factorization {
        self
    }
    fn one(self: &Rc<Self>) -> FpNumber {
        (1, Rc::clone(self))
    }
}

impl SemigroupElem for FpNumber {
    type Group = FpStar;
    fn is_one(&self) -> bool {
        self.0 == 1
    }

    fn group(&self) -> &Rc<FpStar> {
        &self.1
    }

    fn multiply(&mut self, other: &FpNumber) {
        self.0 = long_multiply(self.0, other.0, self.1.p());
    }

    fn square(&mut self) {
        self.0 = long_multiply(self.0, self.0, self.1.p());
    }
}

impl SylowDecomposable for FpStar {
    fn find_sylow_generator(self: &Rc<Self>, d: &(u128, u128)) -> FpNumber {
        match d {
            (2,1) => (self.size().value, Rc::clone(self)),
            _ => {
                (1..self.size().value)
                    .map(|i| (standard_affine_shift(self.p(), i), Rc::clone(self)))
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
        let p = Rc::new(Factorization {
            value: 6,
            prime_powers: vec![(2,1), (3,1)]
        });
        let one = p.one();
        assert!(one.is_one());
    }

    #[test]
    fn multiplies() {
        let p = Rc::new(Factorization {
            value: 6,
            prime_powers: vec![(2,1), (3,1)]
        });
        let mut x = (3, Rc::clone(&p));
        let five = (5, Rc::clone(&p));
        x.multiply(&five);
        assert_eq!(1, x.0);
        assert_eq!(5, five.0);
    }

    #[test]
    fn squares() {
        let p = Rc::new(Factorization {
            value: 6,
            prime_powers: vec![(2,1), (3,1)]
        });
        let mut x = (3, Rc::clone(&p));
        x.square();
        assert_eq!(2, x.0);

        let mut x = (2, Rc::clone(&p));
        x.square();
        assert_eq!(4, x.0);
    }

    #[test]
    fn powers_up() {
        let p = Rc::new(Factorization {
            value: 6,
            prime_powers: vec![(2,1), (3,1)]
        });
        let mut x = (2, Rc::clone(&p));
        x.pow(5);
        assert_eq!(4, x.0);

        let mut x = (3, Rc::clone(&p));
        x.pow(3);
        assert_eq!(6, x.0);

        let mut x = (5, Rc::clone(&p));
        x.pow(p.size().value);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let fp = Rc::new(Factorization {
            value: 1_000_000_000_000_000_124_398,
            prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        });
        let mut x = (3, Rc::clone(&fp));
        x.pow(fp.size().value);
        assert!(x.is_one());
    }

    #[test]
    fn sylow_one_is_one() {
        let fp = Rc::new(Factorization {
            value: 12,
            prime_powers: vec![(2,2), (3,1)]
        });
        let g = Rc::new(SylowDecomp::new(&fp));
        let one = g.one();
        assert!(one.is_one());
    }

    #[test]
    fn sylow_finds_generators() {
        let fp = Rc::new(Factorization {
            value: 19,
            prime_powers: vec![(2, 2), (5, 1)]
        });
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().factor(i);
            test_is_generator_small::<FpStar>(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp = Rc::new(Factorization {
            value: 1_000_000_000_000_000_124_398,
            prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        });
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().prime_powers[i];
            test_is_generator_big::<FpStar>(gen, d);
        }
    }

    #[test]
    fn sylow_order() {
        let fp = Rc::new(Factorization {
            value: 12,
            prime_powers: vec![(2,2), (3,1)]
        });
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 1..13 {
            let mut x = SylowElem {
                group: Rc::clone(&g),
                coords: vec![i % 4, i % 3]
            };
            x.pow(x.order().value);
            assert!(x.is_one());
        }
    }

    #[test]
    fn sylow_order_big() {
        let fp = Rc::new(Factorization {
            value: 1_000_000_000_000_000_124_398,
            prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        });
        let g = Rc::new(SylowDecomp::new(&fp));
        let n = 123456789;
        let mut x = SylowElem {
            group: Rc::clone(&g),
            coords: g.size().prime_powers.iter()
                .map(|(p,d)| n % intpow(*p, *d, 0))
                .collect()
        };
        let or = x.order().value;
        x.pow(or);
        assert!(x.is_one());
    }
}

