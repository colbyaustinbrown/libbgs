use either::*;

use crate::numbers::factorization::*;
use crate::numbers::sylow::*;
use crate::numbers::fp::*;

#[derive(PartialEq, Eq, Debug)]
pub struct QuadField<const P: u128> {
    r: u128
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

impl<const P: u128> QuadField<P> {
    pub fn new(r: u128) -> QuadField<P> {
        QuadField { r }
    }

    pub fn make() -> QuadField<P> {
        QuadField::new(find_nonresidue(P))
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
        r = r.invert(&fp);
        x = x.multiply(&r, &fp);
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

    fn is_one(&self, _f: &QuadField<P>) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&self, other: &QuadNum<P>, f: &QuadField<P>) -> QuadNum<P> {
        QuadNum {
            a0: (long_multiply(self.a0, other.a0, P) + long_multiply(self.a1, long_multiply(other.a1, f.r(), P), P)) % P,
            a1: (long_multiply(self.a1, other.a0, P) + long_multiply(self.a0, other.a1, P)) % P
        }
    }

    fn square(&self, f: &QuadField<P>) -> QuadNum<P> {
        QuadNum {
            a0: (long_multiply(self.a0, self.a0, P) + long_multiply(self.a1, long_multiply(self.a1, f.r(), P), P)) % P,
            a1: (long_multiply(self.a1, self.a0, P) + long_multiply(self.a0, self.a1, P)) % P
        }
    }
}

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
                let p = self.steinitz(j);
                p.pow(pow, self)
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
        let f49 = QuadField::<7>::make();
        let one = f49.one();
        assert!(one.is_one(&f49));
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let f49 = QuadField::<7>::make();
        for i in 2..7 {
            assert_ne!((i * i) % 7, f49.r());
        }
    }

    #[test]
    fn calculates_r_big() {
        QuadField::<BIG_P>::make();
    }

    #[test]
    fn powers_up() {
        let f49 = QuadField::<7>::make();
        let mut x = QuadNum::from_ints(3, 4);
        x = x.pow(48, &f49);
        assert!(x.is_one(&f49));
    }

    #[test]
    fn powers_up_big() {
        let fp2 = QuadField::<BIG_P>::make();
        let mut x = QuadNum {
            a0: 3,
            a1: 5
        };
        x = x.pow(BIG_P - 1, &fp2);
        x = x.pow(BIG_P + 1, &fp2);
        println!("{x:?}");
        assert!(x.is_one(&fp2));
    }

    #[test]
    fn finds_sqrt() {
        let fp2 = QuadField::<BIG_P>::make();
        for i in 3..1003 {
            let mut x = fp2.int_sqrt(i);
            let y = x.clone();
            assert_ne!(x, i);
            x = x.multiply(&y, &fp2);
            assert_eq!(x, i);
        }
    }

    #[test]
    fn sylow_finds_generators() {
        let pplusone = Factorization::new(
            vec![(2, 1), (3, 2)]
        );
        let f289 = QuadField::<17>::make();
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
        let fp2 = QuadField::<BIG_P>::make();
        let g = SylowDecomp::new(&fp2, fp2_fact);
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.factors()[i];
            test_is_generator_big(gen, d, &fp2);
        }
    }
}

