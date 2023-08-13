use either::*;

use crate::numbers::factorization::*;
use crate::numbers::sylow::*;
use crate::numbers::fp::*;

#[derive(PartialEq, Eq, Debug)]
pub struct QuadField<const P: u128> {
    pplusone: Factorization,
    r: u128
}

impl<const P: u128> QuadField<P> {
    pub fn pplusone(&self) -> &Factorization {
        &self.pplusone
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct QuadNum<const P: u128> {
    pub a0: u128,
    pub a1: u128
}

impl<const P: u128> Semigroup for QuadField<P> {
    type Elem = QuadNum<P>;

    fn size(&self) -> u128 {
        P + 1
    }
    fn one(&self) -> QuadNum<P> {
        QuadNum {
            a0: 1,
            a1: 0
        }
    }
}

impl<const P: u128> Factored for QuadField<P> {
    fn factors(&self) -> &Factorization {
        &self.pplusone
    }
}

impl<const P: u128> QuadField<P> {
    pub fn new(pminusone: FpStar<P>, pplusone: Factorization, r: u128) -> QuadField<P> {
        QuadField {
            pplusone,
            r
        }
    }

    pub fn make(pminusone: FpStar<P>, pplusone: Factorization) -> QuadField<P> {
        QuadField::new(pminusone, pplusone, find_nonresidue(P))
    }

    pub fn steinitz(&self, i: u128) -> QuadNum<P> {
        QuadNum::from_ints(i % P, i / P)
    }

    pub fn int_sqrt_either(&self, x: u128) -> Either<QuadNum<P>, FpNum<P>> {
        let fp = FpStar::<P> {};
        let mut x = fp.from_int(x);
        if let Some(y) = x.int_sqrt() {
            return Right(fp.from_int(y.value));
        }

        let mut r = fp.from_int(self.r);
        r.invert(&fp);
        x.multiply(&r, &fp);
        let a1 = x.int_sqrt().unwrap();
        Left(QuadNum {
            a0: 0,
            a1: a1.value
        })
    }

    pub fn int_sqrt(&self, x: u128) -> QuadNum<P> {
        self.int_sqrt_either(x).left_or_else(|n| QuadNum::from_ints(n.value, 0))
    }

    pub fn r(&self) -> u128 {
        self.r
    }

    pub fn change_r(&self, r: u128) -> QuadField<P> {
        QuadField {
            r,
            pplusone: self.pplusone.clone()
        }
    }
}

impl<const P: u128> QuadNum<P> {
    pub fn is_zero(&self) -> bool {
        self.a0 == 0 && self.a1 == 0
    }

    pub fn from_ints(a0: u128, a1: u128) -> Self {
        QuadNum { a0, a1 }
    }

    pub fn add(&mut self, other: QuadNum<P>) {
        self.a0 = (self.a0 + other.a0) % P;
        self.a1 = (self.a1 + other.a1) % P;
    }
}

impl<const P: u128> SemigroupElem for QuadNum<P> {
    type Group = QuadField<P>;

    fn is_one(&self) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&mut self, other: &QuadNum<P>, f: &QuadField<P>) {
        let a0_old = self.a0;
        self.a0 = (long_multiply(self.a0, other.a0, P) + long_multiply(self.a1, long_multiply(other.a1, f.r(), P), P)) % P;
        self.a1 = (long_multiply(self.a1, other.a0, P) + long_multiply(a0_old, other.a1, P)) % P;
    }

    fn square(&mut self, f: &QuadField<P>) {
        let a0_old = self.a0;
        self.a0 = (long_multiply(self.a0, self.a0, P) + long_multiply(self.a1, long_multiply(self.a1, f.r(), P), P)) % P;
        self.a1 = (long_multiply(self.a1, a0_old, P) + long_multiply(a0_old, self.a1, P)) % P;
    }
}
impl<const P: u128> FactoredElem<QuadField<P>> for QuadNum<P> {}

impl<const P: u128> PartialEq<u128> for QuadNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.a0 == *other && self.a1 == 0
    }
}

impl<const P: u128> SylowDecomposable for QuadField<P> {
    fn find_sylow_generator(&self, i: usize, fact: &Factorization) -> QuadNum<P> {
        let pow = P - 1;
        // should be self.p * self.p, but maybe this works?
        (1..P * 2)
            .map(|i| {
                let j = standard_affine_shift(P * 2, i);
                let mut p = self.steinitz(j);
                p.pow(pow, self);
                p
            })
            .find_map(|c| self.is_sylow_generator(&c, fact[i]))
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    fn one_is_one() {
        let f49 = QuadField::make(
            FpStar::<7> {},
            Factorization ::new(vec![(2,3)])
        );
        let one = f49.one();
        assert!(one.is_one());
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let f49 = QuadField::make(
            FpStar::<7> {},
            Factorization ::new(vec![(2,3)])
        );
        for i in 2..7 {
            assert_ne!((i * i) % 7, f49.r());
        }
    }

    #[test]
    fn powers_up() {
        let f49 = QuadField::make(
            FpStar::<7> {},
            Factorization ::new(vec![(2,3)])
        );
        let mut x = QuadNum::from_ints(3, 4);
        x.pow(48, &f49);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let fp2 = QuadField::make(
            FpStar::<BIG_P> {},
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
        let mut x = QuadNum {
            a0: 3,
            a1: 5
        };
        x.pow(BIG_P - 1, &fp2);
        x.pow(BIG_P + 1, &fp2);
        println!("{x:?}");
        assert!(x.is_one());
    }

    #[test]
    fn finds_sqrt() {
        let fp2 = QuadField::make(
            FpStar::<BIG_P> {},
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
        let f289 = QuadField::make(
            FpStar::<17> {},
            pplusone.clone()
        );
        let g = SylowDecomp::new(&f289, pplusone);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors().factor(i);
            test_is_generator_small(gen, d, &f289);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fp2_fact = Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]);
        let fp2 = QuadField::make(
            FpStar::<BIG_P> {},
            fp2_fact.clone()
        );
        let g = SylowDecomp::new(&fp2, fp2_fact);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big(gen, d, &fp2);
        }
    }
}

