use std::marker::PhantomData;

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

    fn from_ints(&self, a0: u128, a1: u128) -> QuadNum<Self> {
        QuadNum {
            _f: PhantomData,
            a0,
            a1
        }
    }

    fn steinitz(&self, i: u128) -> QuadNum<Self> {
        Self::from_ints(self, i % self.p(), i / self.p())
    }
}


#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldExt<'a> {
    pminusone: &'a FpStar,
    pplusone: Factorization,
    r: u128
}

impl<'a> QuadFieldExt<'a> {
    pub fn pminusone(&self) -> &'a FpStar { 
        self.pminusone 
    }
    pub fn pplusone(&self) -> &Factorization {
        &self.pplusone
    }
}

#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldBare {
    p: u128,
    r: u128
}

#[derive(PartialEq, Eq, Debug)]
pub struct QuadNum<F: QuadField> {
    _f: PhantomData<F>,
    pub a0: u128,
    pub a1: u128
}

pub type QuadNumExt<'a> = QuadNum<QuadFieldExt<'a>>;

impl QuadFieldBare {
    pub fn new(p: u128, r: u128) -> QuadFieldBare {
        QuadFieldBare {
            p,
            r
        }
    }
}

impl<'a> Semigroup for QuadFieldExt<'a> {
    type Elem = QuadNum<QuadFieldExt<'a>>;

    fn size(&self) -> u128 {
        self.pplusone.value()
    }
    fn one(&self) -> QuadNum<QuadFieldExt<'a>> {
        QuadNum {
            _f: PhantomData,
            a0: 1,
            a1: 0
        }
    }
}

impl Semigroup for QuadFieldBare {
    type Elem = QuadNum<QuadFieldBare>;

    fn size(&self) -> u128 {
        self.p - 1
    }
    fn one(&self) -> QuadNum<QuadFieldBare> {
        QuadNum {
            _f: PhantomData,
            a0: 1,
            a1: 0
        }
    }
}

impl<'a> Factored for QuadFieldExt<'a> {
    fn factors(&self) -> &Factorization {
        &self.pplusone
    }
}

impl<'a> QuadFieldExt<'a> {
    pub fn new(pminusone: &'a FpStar, pplusone: Factorization, r: u128) -> QuadFieldExt {
        QuadFieldExt {
            pminusone: pminusone,
            pplusone,
            r
        }
    }

    pub fn make(pminusone: &FpStar, pplusone: Factorization) -> QuadFieldExt {
        let p = pplusone.value() - 1;
        QuadFieldExt::new(pminusone, pplusone, find_nonresidue(p))
    }

    pub fn int_sqrt_either(&self, x: u128) -> Either<QuadNum<QuadFieldExt>, FpNum> {
        let mut x = self.pminusone.from_int(x);
        let fp = self.pminusone;
        if let Some(y) = x.int_sqrt(fp) {
            return Right(self.pminusone.from_int(y.value()));
        }

        let mut r = self.pminusone.from_int(self.r);
        r.invert(fp);
        x.multiply(&r, fp);
        let a1 = x.int_sqrt(fp).unwrap();
        Left(QuadNum {
            _f: PhantomData,
            a0: 0,
            a1: a1.value()
        })
    }

    pub fn int_sqrt(&self, x: u128) -> QuadNum<QuadFieldExt> {
        self.int_sqrt_either(x).left_or_else(|n| self.from_ints(n.value(), 0))
    }
}

impl<'a> QuadField for QuadFieldExt<'a> {
    fn r(&self) -> u128 {
        self.r
    }

    fn change_r(&self, r: u128) -> QuadFieldExt<'a> {
        QuadFieldExt {
            r: r,
            pminusone: &self.pminusone,
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

impl<F: QuadField> Clone for QuadNum<F> {
    fn clone(&self) -> QuadNum<F> {
        QuadNum {
            _f: PhantomData,
            a0: self.a0,
            a1: self.a1
        }
    }
}

impl<F> QuadNum<F> 
where F: QuadField<Elem = Self> {
    fn is_zero(&self) -> bool {
        self.a0 == 0 && self.a1 == 0
    }

    fn from_int(a0: u128) -> QuadNum<F> {
        QuadNum {
            _f: PhantomData,
            a0,
            a1: 0
        }
    }
}

impl<F> SemigroupElem for QuadNum<F> 
where F: QuadField<Elem = Self> {
    type Group = F;

    fn is_one(&self, _: &F) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&mut self, other: &QuadNum<F>, f: &F) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, other.a0, p) + long_multiply(a0_old, other.a1, p)) % p;
    }

    fn square(&mut self, f: &F) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, a0_old, p) + long_multiply(a0_old, self.a1, p)) % p;
    }
}
impl<'a> FactoredElem<QuadFieldExt<'a>> for QuadNum<QuadFieldExt<'a>> {}

impl<F: QuadField> PartialEq<u128> for QuadNum<F> {
    fn eq(&self, other: &u128) -> bool {
        self.a0 == *other && self.a1 == 0
    }
}

impl<'a> SylowDecomposable for QuadFieldExt<'a> {
    fn find_sylow_generator(&self, i: usize) -> QuadNum<QuadFieldExt<'a>> {
        let pow = self.pminusone.value();
        // should be self.p * self.p, but maybe this works?
        (1..self.p() * 2)
            .map(|i| {
                let j = standard_affine_shift(self.p() * 2, i);
                let mut p = self.steinitz(j);
                p.pow(pow, &self);
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

    #[test]
    fn one_is_one() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadFieldExt::make(
            &p7,
            Factorization ::new(vec![(2,3)])
        );
        let one = f49.one();
        assert!(one.is_one(&f49));
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadFieldExt::make(
            &p7,
            Factorization ::new(vec![(2,3)])
        );
        for i in 2..f49.p() {
            assert_ne!((i * i) % f49.p(), f49.r());
        }
    }

    #[test]
    fn powers_up() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadFieldExt::make(
            &p7,
            Factorization ::new(vec![(2,3)])
        );
        let mut x = f49.from_ints(3, 4);
        x.pow(48, &f49);
        assert!(x.is_one(&f49));
    }

    #[test]
    fn powers_up_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadFieldExt::make(
            &fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        let mut x = QuadNum {
            _f: PhantomData,
            a0: 3,
            a1: 5
        };
        x.pow(fp2.pminusone.value(), &fp2);
        x.pow(fp2.pplusone.value(), &fp2);
        assert!(x.is_one(&fp2));
    }

    #[test]
    fn finds_sqrt() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadFieldExt::make(
            &fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        for i in 3..1003 {
            let mut x = fp2.int_sqrt(i);
            let y = x.clone();
            assert_ne!(x, i);
            x.multiply(&y, &fp2);
            assert_eq!(x, i);
        }
    }

    #[test]
    fn sylow_finds_generators() {
        let pplusone = Factorization::new(
            vec![(2, 1), (3, 2)]
        );
        let f17 = Factorization::new(vec![(2,4)]);
        let f289 = QuadFieldExt::make(
            &f17,
            pplusone.clone()
        );
        let g = SylowDecomp::new(&f289);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors().factor(i);
            test_is_generator_small::<QuadFieldExt>(gen, d, &f289);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadFieldExt::make(
            &fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        let g = SylowDecomp::new(&fp2);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big::<QuadFieldExt>(gen, d, &fp2);
        }
    }
}

