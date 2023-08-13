pub use crate::numbers::semigroup::*;
use crate::util::long_multiply;
use crate::numbers::sylow::*;
use crate::util::*;
use crate::numbers::factorization::*;
use crate::numbers::quad_field_small::*;
pub use crate::numbers::group::*;

#[derive(PartialEq, Clone, Copy, Debug, Eq)]
pub struct FpStar<const P: u128> {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FpNum<const P: u128> {
    pub value: u128
}

impl<const P: u128> FpStar<P> {
    pub fn from_int(&self, value: u128) -> FpNum<P> {
        FpNum {
            value
        }
    }

    pub fn p(&self) -> u128 {
        P
    }
}

impl<const P: u128> FpNum<P> {
    pub fn add(&mut self, other: &FpNum<P>) {
        self.value += other.value;
        self.value %= P;
    }

    pub fn sub(&mut self, other: &FpNum<P>) {
        while self.value < other.value {
            self.value += P;
        }
        self.value -= other.value;
    }

    pub fn int_sqrt(&self) -> Option<FpNum<P>> {
        let fp = FpStar::<P> {};
        if self.value == 0 { 
            return Some(fp.from_int(0));
        }

        // Cipolla's algorithm
        let l = legendre(self.value, P);
        if l == 0 { 
            return Some(fp.from_int(0));
        } else if l == P - 1{
            return None;
        }
        let mut i = 1;
        let (a,r) = loop {
            let a = standard_affine_shift(P, i);
            let r = (intpow(a, 2, P) + P - self.value) % P;
            if legendre(r, P) == (P - 1) { break (a,r); }
            i += 1;
        };
        let fp2 = QuadFieldSml::new(P, r);
        let mut x = QuadNumSml::from_ints(a, 1);
        x.pow((P + 1) / 2, &fp2);
        Some(fp.from_int(x.a0))
    }
}

impl<const P: u128> PartialEq<u128> for FpNum<P> {
    fn eq(&self, other: &u128) -> bool {
        self.value == *other
    }
}

impl<const P: u128> Semigroup for FpStar<P> {
    type Elem = FpNum<P>;

    fn size(&self) -> u128 {
        P - 1
    }

    fn one(&self) -> FpNum<P> {
        self.from_int(1)
    }
}

impl<const P: u128> Group for FpStar<P> {}

impl<const P: u128> SylowDecomposable for FpStar<P> {
    fn find_sylow_generator(&self, i: usize, fact: &Factorization) -> FpNum<P> {
        match fact[i] {
            (2,1) => self.from_int(self.size()),
            _ => (1..self.size())
                .map(|i| self.from_int(standard_affine_shift(P, i)))
                .inspect(|i| println!("i is {i:?}"))
                .find_map(|c| self.is_sylow_generator(&c, fact[i]))
                .unwrap()
        }
    }
}

impl<const P: u128> SemigroupElem for FpNum<P> {
    type Group = FpStar<P>;
    fn is_one(&self) -> bool {
        self.value == 1
    }

    fn multiply(&mut self, other: &FpNum<P>, _: &FpStar<P>) {
        self.value = long_multiply(self.value, other.value, P);
    }

    fn square(&mut self, _: &FpStar<P>) {
        self.value = long_multiply(self.value, self.value, P);
    }
}

impl<const P: u128> GroupElem for FpNum<P> {
    fn invert(&mut self, g: &FpStar<P>) {
        self.pow(P - 2, g); 
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::sylow::tests::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    fn one_is_one() {
        let p = FpStar::<7> {};
        let one = p.one();
        assert!(one.is_one());
    }

    #[test]
    fn multiplies() {
        let p = FpStar::<7> {};
        let mut x = p.from_int(3);
        assert_eq!(3, x.value);
        let five = p.from_int(5);
        x.multiply(&five, &p);
        assert_eq!(1, x.value);
        assert_eq!(5, five.value);
    }

    #[test]
    fn squares() {
        let p = FpStar::<7> {};
        let mut x = p.from_int(3);
        x.square(&p);
        assert_eq!(2, x.value);

        let mut x = p.from_int(2);
        x.square(&p);
        assert_eq!(4, x.value);
    }

    #[test]
    fn powers_up() {
        let p = FpStar::<7> {};
        let mut x = p.from_int(2);
        x.pow(5, &p);
        assert_eq!(4, x.value);

        let mut x = p.from_int(3);
        x.pow(3, &p);
        assert_eq!(6, x.value);

        let mut x = p.from_int(5);
        x.pow(p.p() - 1, &p);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let p = FpStar::<BIG_P> {};
        let mut x = p.from_int(3);
        x.pow(BIG_P - 1, &FpStar::<BIG_P> {});
        println!("x is {x:?}");
        assert!(x.is_one());
    }

    #[test]
    fn sylow_one_is_one() {
        let p = Factorization::new(vec![(2,2), (3,1)]);
        let g = SylowDecomp::new(&FpStar::<13> {}, p.clone());
        let one = g.one();
        assert!(one.is_one());
    }

    #[test]
    fn sylow_finds_generators() {
        let p = Factorization::new(vec![(2, 2), (7, 1)]);
        let g = SylowDecomp::new(&FpStar::<29> {}, p.clone());
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.fact.factor(i);
            test_is_generator_small::<FpStar<29>>(gen, d, &FpStar::<29> {});
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let fact = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&FpStar::<BIG_P> {}, fact.clone());
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = fact[i];
            test_is_generator_big::<FpStar<BIG_P>>(gen, d, &FpStar::<BIG_P> {});
        }
    }

    #[test]
    fn sylow_order() {
        let fact = Factorization::new(vec![(2,2), (3,1)]);
        let g = SylowDecomp::new(&FpStar::<13> {}, fact.clone());
        for i in 1..13 {
            let mut x = SylowElem::new(vec![i % 4, i % 3]);
            x.pow(x.order(&g).value(), &g);
            assert!(x.is_one());
        }
    }

    #[test]
    fn sylow_order_big() {
        let fact = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&FpStar::<BIG_P> {}, fact);
        let n = 123456789;
        let mut x = SylowElem::new(
            g.factors().as_array().iter()
                .map(|(p,d)| n % intpow(*p, *d, 0))
                .collect()
        );
        let or = x.order(&g).value();
        x.pow(or, &g);
        assert!(x.is_one());
    }

    #[test]
    fn calculates_square_roots() {
        let p = FpStar::<13> {};
        let mut nonresidues = 0;
        for x in (1..13).map(|i| p.from_int(i)) {
            match x.int_sqrt() {
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
        let p = FpStar::<13> {};
        for i in 2..13 {
            let mut x = p.from_int(i);
            let y = x.clone();
            x.invert(&p);
            assert!(!x.is_one());
            x.multiply(&y, &p);
            assert!(x.is_one());
        }
    }
}

