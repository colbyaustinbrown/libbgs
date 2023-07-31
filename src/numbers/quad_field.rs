use either::*;

use crate::numbers::factorization::*;
use crate::numbers::sylow::*;
use crate::numbers::fp::*;
use crate::numbers::group::*;

#[derive(PartialEq, Eq, Debug)]
pub struct QuadField {
    pminusone: FpStar,
    pplusone: Factorization,
    r: u128
}

impl QuadField {
    pub fn pminusone(&self) -> &FpStar { 
        &self.pminusone 
    }
    pub fn pplusone(&self) -> &Factorization {
        &self.pplusone
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct QuadNum {
    pub a0: u128,
    pub a1: u128
}

impl Semigroup for QuadField {
    type Elem = QuadNum;

    fn size(&self) -> u128 {
        self.pplusone.value()
    }
    fn one(&self) -> QuadNum {
        QuadNum {
            a0: 1,
            a1: 0
        }
    }
}

impl Factored for QuadField {
    fn factors(&self) -> &Factorization {
        &self.pplusone
    }
}

impl QuadField {
    pub fn new(pminusone: FpStar, pplusone: Factorization, r: u128) -> QuadField {
        QuadField {
            pminusone: pminusone,
            pplusone,
            r
        }
    }

    pub fn make(pminusone: FpStar, pplusone: Factorization) -> QuadField {
        let p = pplusone.value() - 1;
        QuadField::new(pminusone, pplusone, find_nonresidue(p))
    }

    pub fn steinitz(&self, i: u128) -> QuadNum {
        QuadNum::from_ints(i % self.p(), i / self.p())
    }

    pub fn int_sqrt_either(&self, x: u128) -> Either<QuadNum, FpNum> {
        let mut x = self.pminusone.from_int(x);
        let fp = &self.pminusone;
        if let Some(y) = x.int_sqrt(fp) {
            return Right(self.pminusone.from_int(y.value()));
        }

        let mut r = self.pminusone.from_int(self.r);
        r.invert(fp);
        x.multiply(&r, fp);
        let a1 = x.int_sqrt(fp).unwrap();
        Left(QuadNum {
            a0: 0,
            a1: a1.value()
        })
    }

    pub fn int_sqrt(&self, x: u128) -> QuadNum {
        self.int_sqrt_either(x).left_or_else(|n| QuadNum::from_ints(n.value(), 0))
    }

    pub fn r(&self) -> u128 {
        self.r
    }

    pub fn change_r(&self, r: u128) -> QuadField {
        QuadField {
            r: r,
            pminusone: self.pminusone.clone(),
            pplusone: self.pplusone.clone()
        }
    }

    pub fn p(&self) -> u128 {
        self.pminusone.value() + 1
    }
}

impl QuadNum {
    pub fn is_zero(&self) -> bool {
        self.a0 == 0 && self.a1 == 0
    }

    pub fn from_ints(a0: u128, a1: u128) -> Self {
        QuadNum { a0, a1 }
    }
}

impl SemigroupElem for QuadNum {
    type Group = QuadField;

    fn is_one(&self, _f: &QuadField) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&mut self, other: &QuadNum, f: &QuadField) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, other.a0, p) + long_multiply(a0_old, other.a1, p)) % p;
    }

    fn square(&mut self, f: &QuadField) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, a0_old, p) + long_multiply(a0_old, self.a1, p)) % p;
    }
}
impl FactoredElem<QuadField> for QuadNum {}

impl PartialEq<u128> for QuadNum {
    fn eq(&self, other: &u128) -> bool {
        self.a0 == *other && self.a1 == 0
    }
}

impl SylowDecomposable for QuadField {
    fn find_sylow_generator(&self, i: usize) -> QuadNum {
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
    use crate::numbers::sylow::tests::*;

    #[test]
    fn one_is_one() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadField::make(
            p7,
            Factorization ::new(vec![(2,3)])
        );
        let one = f49.one();
        assert!(one.is_one(&f49));
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadField::make(
            p7,
            Factorization ::new(vec![(2,3)])
        );
        for i in 2..f49.p() {
            assert_ne!((i * i) % f49.p(), f49.r());
        }
    }

    #[test]
    fn powers_up() {
        let p7 = Factorization::new(vec![(2,1), (3,1)]);
        let f49 = QuadField::make(
            p7,
            Factorization ::new(vec![(2,3)])
        );
        let mut x = QuadNum::from_ints(3, 4);
        x.pow(48, &f49);
        assert!(x.is_one(&f49));
    }

    #[test]
    fn powers_up_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadField::make(
            fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        let mut x = QuadNum {
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
        let fp2 = QuadField::make(
            fp,
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
        let f289 = QuadField::make(
            f17,
            pplusone.clone()
        );
        let g = SylowDecomp::new(&f289);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors().factor(i);
            test_is_generator_small::<QuadField>(gen, d, &f289);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadField::make(
            fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        let g = SylowDecomp::new(&fp2);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big::<QuadField>(gen, d, &fp2);
        }
    }
}

