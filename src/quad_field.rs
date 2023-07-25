use std::rc::Rc;

use either::*;

use crate::factorization::*;
use crate::sylow::*;
use crate::fp::*;
use crate::group::*;

pub trait QuadField: Semigroup + Sized + PartialEq + Eq + std::fmt::Debug 
where Self: Semigroup {
    fn p(&self) -> u128;
    fn r(&self) -> u128;
    fn change_r(&self, r: u128) -> Self;

    fn from_ints(self: &Rc<Self>, a0: u128, a1: u128) -> QuadNumber<Self> {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0,
            a1
        }
    }

    fn steinitz(self: &Rc<Self>, i: u128) -> QuadNumber<Self> {
        Self::from_ints(self, i % self.p(), i / self.p())
    }
}


#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldExt {
    pminusone: Rc<FpStar>,
    pplusone: Factorization,
    r: u128
}

#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldBare {
    p: u128,
    r: u128
}

impl QuadFieldBare {
    pub fn new(p: u128, r: u128) -> QuadFieldBare {
        QuadFieldBare {
            p,
            r
        }
    }
}

impl Semigroup for QuadFieldExt {
    type Elem = QuadNumber;

    fn size(&self) -> u128 {
        self.pplusone.value()
    }
    fn one(self: &Rc<Self>) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0: 1,
            a1: 0
        }
    }
}

impl Semigroup for QuadFieldBare {
    type Elem = QuadNumber<QuadFieldBare>;

    fn size(&self) -> u128 {
        self.p - 1
    }
    fn one(self: &Rc<Self>) -> QuadNumber<QuadFieldBare> {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0: 1,
            a1: 0
        }
    }
}

impl Factored for QuadFieldExt {
    fn factors(&self) -> &Factorization {
        &self.pplusone
    }
}

impl QuadFieldExt {
    pub fn new(pminusone: &Rc<FpStar>, pplusone: Factorization, r: u128) -> QuadFieldExt {
        QuadFieldExt {
            pminusone: Rc::clone(pminusone),
            pplusone,
            r
        }
    }

    pub fn make(pminusone: &Rc<FpStar>, pplusone: Factorization) -> QuadFieldExt {
        let p = pplusone.value() - 1;
        QuadFieldExt::new(pminusone, pplusone, find_nonresidue(p))
    }

    pub fn int_sqrt_either(self: &Rc<Self>, x: u128) -> Either<QuadNumber, FpNumber> {
        let mut x = self.pminusone.from_int(x);
        if let Some(y) = x.int_sqrt() {
            return Right(self.pminusone.from_int(y.value()));
        }

        let mut r = self.pminusone.from_int(self.r);
        r.invert();
        x.multiply(&r);
        let a1 = x.int_sqrt().unwrap();
        Left(QuadNumber {
            subgroup: Rc::clone(self),
            a0: 0,
            a1: a1.value()
        })
    }

    pub fn int_sqrt(self: &Rc<Self>, x: u128) -> QuadNumber {
        self.int_sqrt_either(x).left_or_else(|n| self.from_ints(n.value(), 0))
    }
}

impl QuadField for QuadFieldExt {
    fn r(&self) -> u128 {
        self.r
    }

    fn change_r(&self, r: u128) -> QuadFieldExt {
        QuadFieldExt {
            r: r,
            pminusone: Rc::clone(&self.pminusone),
            pplusone: self.pplusone.clone()
        }
    }

    fn p(&self) -> u128 {
        self.pminusone.value() + 1
    }
}

impl QuadField for QuadFieldBare {
    fn r(&self) -> u128 {
        self.r
    }

    fn p(&self) -> u128 {
        self.p
    }

    fn change_r(&self, r: u128) -> QuadFieldBare {
        QuadFieldBare {
            p: self.p,
            r
        }
    }
}

#[derive(PartialEq, Eq, Debug)]
pub struct QuadNumber<F: QuadField = QuadFieldExt> {
    subgroup: Rc<F>,
    pub a0: u128,
    pub a1: u128
}

impl<F: QuadField> Clone for QuadNumber<F> {
    fn clone(&self) -> QuadNumber<F> {
        QuadNumber {
            subgroup: Rc::clone(&self.subgroup),
            a0: self.a0,
            a1: self.a1
        }
    }
}

impl<F: QuadField> QuadNumber<F> {
    fn is_zero(&self) -> bool {
        self.a0 == 0 && self.a1 == 0
    }

    fn from_int(subgroup: Rc<F>, a0: u128) -> QuadNumber<F> {
        QuadNumber {
            subgroup,
            a0,
            a1: 0
        }
    }
}

impl<F> SemigroupElem for QuadNumber<F> 
where F: QuadField<Elem = Self> {
    type Group = F;

    fn is_one(&self) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn group(&self) -> &Rc<F> {
        &self.subgroup
    }

    fn multiply(&mut self, other: &QuadNumber<F>) {
        let f = &self.subgroup;
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, other.a0, p) + long_multiply(a0_old, other.a1, p)) % p;
    }

    fn square(&mut self) {
        let f = &self.subgroup;
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, a0_old, p) + long_multiply(a0_old, self.a1, p)) % p;
    }
}
impl FactoredElem<QuadFieldExt> for QuadNumber {}

impl<F: QuadField> PartialEq<u128> for QuadNumber<F> {
    fn eq(&self, other: &u128) -> bool {
        self.a0 == *other && self.a1 == 0
    }
}

impl SylowDecomposable for QuadFieldExt {
    fn find_sylow_generator(self: &Rc<Self>, i: usize) -> QuadNumber<QuadFieldExt> {
        let pow = self.pminusone.value();
        // should be self.p * self.p, but maybe this works?
        (1..self.p() * 2)
            .map(|i| {
                let j = standard_affine_shift(self.p() * 2, i);
                let mut p = self.steinitz(j);
                p.pow(pow);
                p
            })
            .find_map(|c| self.is_sylow_generator(&c, i))
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sylow::tests::*;

    fn p7() -> QuadFieldExt {
        QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2,1), (3,1)]
            )),
            Factorization ::new(
                vec![(2,3)]
            )
        )
    }

    #[test]
    fn one_is_one() {
        let one = Rc::new(p7()).one();
        assert!(one.is_one());
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let f = p7();
        for i in 2..f.p() {
            assert_ne!((i * i) % f.p(), f.r);
        }
    }

    #[test]
    fn powers_up() {
        let mut x = QuadNumber {
            subgroup: Rc::new(p7()),
            a0: 3,
            a1: 4
        };
        x.pow(48);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let fp = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            )),
            Factorization::new(
                vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
            )
        ));
        let mut x = QuadNumber {
            subgroup: Rc::clone(&fp),
            a0: 3,
            a1: 5
        };
        x.pow(fp.pminusone.value());
        x.pow(fp.pplusone.value());
        assert!(x.is_one());
    }

    #[test]
    fn finds_sqrt() {
        let fp = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            )),
            Factorization::new(
                vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
            )
        ));
        for i in 3..1003 {
            let mut x = fp.int_sqrt(i);
            let y = x.clone();
            assert_ne!(x, i);
            x.multiply(&y);
            assert_eq!(x, i);
        }
    }

    #[test]
    fn sylow_finds_generators() {
        let pplusone = Factorization::new(
            vec![(2, 1), (3, 2)]
        );
        let f = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 4)]
            )),
            pplusone.clone()
        ));
        let g = Rc::new(SylowDecomp::new(&f));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors().factor(i);
            test_is_generator_small::<QuadFieldExt>(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let pplusone = Factorization::new(
            vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
        );
        let fp = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            )),
            pplusone.clone()
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big::<QuadFieldExt>(gen, d);
        }
    }
}

