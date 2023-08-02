use crate::numbers::sylow::*;
use crate::util::*;
use crate::numbers::factorization::*;

pub mod flags {
    pub const NONE: u8 = 0b0;
    pub const NO_UPPER_HALF: u8 = 0b1;
}

/*
pub fn sylow_factory<'a, C: SylowDecomposable>(decomp: &'a SylowDecomp<C>, rs: &[u128], mode: Mode) -> impl Iterator<Item = SylowElem<'a, C>> {
    rs.iter()
        .zip(decomp.factors().as_array())
        .enumerate()
        .filter(|(_, (r, _))| **r > 0)
        .map(|(i, (r, (p, d)))| SylowFactory {
            decomp,
            i,
            r: *r,
            stack: vec![(0, intpow(*p, *d - 1, 0), *r)],
            mode
        })
        .multi_cartesian_product()
        .map(|x| {
            let mut res = decomp.one();
            x.iter().for_each(|e| { res.multiply(e, decomp); });
            res
        })
}
*/

struct Gear {
    i: usize,
    coord: u128,
    step: u128,
    r: u128,
    first: bool
}

pub struct SylowFactory<'a, C: SylowDecomposable> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,

    stack: Vec<Gear>
}

impl<'a, C: SylowDecomposable> SylowFactory<'a, C> {
    pub fn new(decomp: &'a SylowDecomp<'a, C>, i: usize, r: u128, mode: u8) -> SylowFactory<'a, C>{
        let (p,d) = decomp.factors()[i];
        SylowFactory {
            decomp,
            mode,
            stack: vec![Gear {
                i,
                coord: 0, 
                step: intpow(p, d - 1, 0), 
                r,
                first: true
            }]
        }
    }
}

impl<'a, C: SylowDecomposable> Iterator for SylowFactory<'a, C> {
    type Item = SylowElem<'a, C>;

    fn next(&mut self) -> Option<SylowElem<'a, C>> {
        loop {
            let l = self.stack.len();
            if l == 0 { return None; }

            let gear = self.stack.remove(l - 1);
            let p = self.decomp.factors()[gear.i].0;

            if gear.r > 0 {
                for j in 1..(p + 1) {
                    if gear.r == 1 && j == p + 1 { continue; }

                    let coord = gear.coord + (p - j) * gear.step;
                    if self.mode & flags::NO_UPPER_HALF != 0 {
                        if coord > self.decomp.factors().factor(gear.i) / 2 {
                                continue;
                        }
                    }

                    self.stack.push(Gear {
                        i: gear.i,
                        coord: gear.coord + (p - j) * gear.step, 
                        step: gear.step / p, 
                        r: gear.r - 1,
                        first: true
                    }); 
                }
            }
            if gear.r == 0 {
                let mut coords = vec![0 ; self.decomp.factors().len()];
                coords[gear.i] = gear.coord;
                return Some(SylowElem::new(coords));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::fp::*;

    #[test]
    pub fn test_make_factory() {
        let fp = Factorization::new(vec![(2, 1), (3, 1)]);
        let g = SylowDecomp::new(&fp);
        let mut factory = SylowFactory::new(&g, 0, 1);
        assert_eq!(factory.next().map(|s| s.to_product(&g).value()), Some(6));
        assert_eq!(factory.next(), None);
    }

    #[test]
    pub fn test_generates_small() {
        let fp = Factorization::new(vec![(2, 2), (3, 1), (5, 1)]);
        let g = SylowDecomp::new(&fp);
        let factory = SylowFactory::new(&g, 0, 1);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 1);
        let mut x = res[0].clone();
        assert!(!x.is_one(&g));
        x.pow(2, &g);
        assert!(x.is_one(&g));

        let factory = SylowFactory::new(&g, 0, 2);
        let mut count = 0;
        for mut x in factory {
            count += 1;
            for _ in 1..3 {
                let y = x.clone();
                assert!(!x.is_one(&g));
                x.multiply(&y, &g);
            }
            assert!(x.is_one(&g));
        }
        assert_eq!(count, 2);

        let factory = SylowFactory::new(&g, 1, 1);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 2);
    }

    #[test]
    pub fn test_generates_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&fp);

        let factory = SylowFactory::new(&g, 3, 2);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 29 * 29 - 29);

        let mut factory = SylowFactory::new(&g, 5, 1);
        let mut x = factory.next();
        assert_eq!(x.as_ref().map(|a| a.is_one(&g)), Some(false));
        x.as_mut().map(|a| {a.pow(705737, &g); a});
        assert_eq!(x.map(|a| a.is_one(&g)), Some(true));
    }
}

