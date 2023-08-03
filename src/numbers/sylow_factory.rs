use crate::numbers::sylow::*;
use crate::util::*;
use crate::numbers::factorization::*;

pub mod flags {
    pub const NONE: u8 = 0x01;
    pub const NO_UPPER_HALF: u8 = 0x02;
}

#[derive(Debug)]
struct Gear<'a, C: SylowDecomposable + std::fmt::Debug> {
    i: usize,
    res: SylowElem<'a, C>,
    step: u128,
    r: u128,
    is_fresh: bool
}

pub struct SylowFactory<'a, C: SylowDecomposable + std::fmt::Debug> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,
    pows: Vec<u128>,
    stack: Vec<Gear<'a, C>>
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> SylowFactory<'a, C> {
    pub fn new(decomp: &'a SylowDecomp<'a, C>, pows: Vec<u128>, mode: u8) -> SylowFactory<'a, C>{
        let first_nonzero = pows.iter()
            .enumerate()
            .find_map(|(i, r)| if *r > 0 { Some(i) } else { None })
            .unwrap_or(0);
        let (p,d) = decomp.factors()[first_nonzero];
        SylowFactory {
            decomp,
            mode,
            stack: vec![Gear {
                i: first_nonzero,
                res: decomp.one(),
                step: intpow(p, d - 1, 0), 
                r: pows[first_nonzero],
                is_fresh: true
            }],
            pows
        }
    }
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> Iterator for SylowFactory<'a, C> {
    type Item = SylowElem<'a, C>;

    fn next(&mut self) -> Option<SylowElem<'a, C>> {
        loop {
            let Some(mut gear) = self.stack.pop() else { return None; };
            let (p,_) = self.decomp.factors()[gear.i];
            let fact = self.decomp.factors().factor(gear.i);

            if gear.r > 0 {
                let start = if gear.r > 1 { 0 } else { 1 };
                for i in start..p { 
                    let mut res = gear.res.clone();
                    res.coords[gear.i] += i * gear.step;

                    if self.mode & flags::NO_UPPER_HALF != 0
                        && gear.is_fresh 
                        && res.coords[gear.i] > fact / 2 { 
                        break; 
                    }

                    self.stack.push(Gear {
                        i: gear.i,
                        res,
                        step: gear.step / p,
                        r: gear.r - 1,
                        is_fresh: gear.is_fresh
                    });
                }
            } else {
                gear.i += 1;

                if gear.i == self.decomp.factors().len() {
                    // println!("gear is {:?}", gear);
                    return Some(gear.res);
                }

                let (p,d) = self.decomp.factors()[gear.i];

                gear.step = intpow(p, d - 1, 0);
                gear.r = self.pows[gear.i];
                gear.is_fresh = false;
                self.stack.push(gear);
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
        let mut factory = SylowFactory::new(&g, vec![1, 0], flags::NONE);
        assert_eq!(factory.next().map(|s| s.to_product(&g).value()), Some(6));
        assert_eq!(factory.next(), None);
    }

    #[test]
    pub fn test_generates_small() {
        let fp = Factorization::new(vec![(2, 2), (3, 1), (5, 1)]);
        let g = SylowDecomp::new(&fp);
        let factory = SylowFactory::new(&g, vec![1, 0, 0], flags::NONE);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 1);
        let mut x = res[0].clone();
        assert!(!x.is_one(&g));
        x.pow(2, &g);
        assert!(x.is_one(&g));

        let factory = SylowFactory::new(&g, vec![2, 0, 0], flags::NONE);
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

        let factory = SylowFactory::new(&g, vec![0, 1, 0], flags::NONE);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 2);
    }

    #[test]
    pub fn test_generates_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&fp);

        let factory = SylowFactory::new(&g, vec![0, 0, 0, 2, 0, 0, 0], flags::NONE);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 29 * 29 - 29);

        let mut factory = SylowFactory::new(&g, vec![0, 0, 0, 0, 0, 1, 0], flags::NONE);
        let mut x = factory.next();
        assert_eq!(x.as_ref().map(|a| a.is_one(&g)), Some(false));
        x.as_mut().map(|a| {a.pow(705737, &g); a});
        assert_eq!(x.map(|a| a.is_one(&g)), Some(true));
    }

    #[test]
    pub fn test_skips_upper_half() {
        let fp = Factorization::new(vec![(2,1), (3,3), (5,1)]);
        let g = SylowDecomp::new(&fp);

        let factory_all = SylowFactory::new(&g, vec![0, 2, 1], flags::NONE);
        let res: Vec<SylowElem<FpStar>> = factory_all.collect();
        println!("{res:?}");
        assert_eq!(res.len(), 24);

        let factory_half = SylowFactory::new(&g, vec![0, 2, 1], flags::NO_UPPER_HALF);
        let res: Vec<SylowElem<FpStar>> = factory_half.collect();
        println!("{res:?}");
        assert_eq!(res.len(), 12);
    }
}

