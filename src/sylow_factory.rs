use std::rc::Rc;

use crate::sylow::*;
use crate::util::*;
use crate::semigroup::*;

#[derive(Debug)]
pub struct SylowFactory<C: SylowDecomposable> {
    pub decomp: Rc<SylowDecomp<C>>,
    pub i: usize,

    stack: Vec<(u128, u128, usize)>
}

impl<C: SylowDecomposable> SylowFactory<C> {
    pub fn new(decomp: &Rc<SylowDecomp<C>>, i: usize, r: usize) -> SylowFactory<C> {
        let (p,d) = decomp.size().prime_powers()[i];
        let step = intpow(p, d - 1, 0);
        SylowFactory {
            decomp: Rc::clone(decomp),
            i,
            stack: vec![(0,step,r)]
        }
    }
}

impl<C: SylowDecomposable> Iterator for SylowFactory<C> {
    type Item = SylowElem<C>;

    fn next(&mut self) -> Option<SylowElem<C>> {
        loop {
            let l = self.stack.len();
            if l == 0 { return None; }

            let p = self.decomp.size().prime_powers()[self.i].0;
            let (x, step, r) = self.stack.remove(l - 1);
            if r == 0 {
                let mut coords = vec![0 ; self.decomp.len()];
                coords[self.i] = x;
                return Some(SylowElem {
                    group: Rc::clone(&self.decomp),
                    coords
                });
            } else {
                for i in 1..p {
                    self.stack.push((x + (p - i) * step, step / p, r - 1)); 
                }
                if r != 1 { self.stack.push((x, step / p, r - 1)); }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factorization::*;
    use crate::fp::*;

    #[test]
    pub fn test_make_factory() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 1), (3, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        let mut factory = SylowFactory::new(&g, 0, 1);
        assert_eq!(factory.next().map(|s| s.to_product().0), Some(6));
        assert_eq!(factory.next(), None);
    }

    #[test]
    pub fn test_generates_small() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 2), (3, 1), (5, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        let factory = SylowFactory::new(&g, 0, 1);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 1);
        let mut x = res[0].clone();
        assert!(!x.is_one());
        x.pow(2);
        assert!(x.is_one());

        let factory = SylowFactory::new(&g, 0, 2);
        let mut count = 0;
        for mut x in factory {
            count += 1;
            for _ in 1..3 {
                let y = x.clone();
                assert!(!x.is_one());
                x.multiply(&y);
            }
            assert!(x.is_one());
        }
        assert_eq!(count, 2);

        let factory = SylowFactory::new(&g, 1, 1);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 2);
    }

    #[test]
    pub fn test_generates_big() {
        let fp = Rc::new(Factorization::new(
            vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        ));
        let g = Rc::new(SylowDecomp::new(&fp));

        let factory = SylowFactory::new(&g, 3, 2);
        let res: Vec<SylowElem<FpStar>> = factory.collect();
        assert_eq!(res.len(), 29 * 29 - 29);

        let mut factory = SylowFactory::new(&g, 5, 1);
        let mut x = factory.next();
        assert_eq!(x.as_ref().map(|a| a.is_one()), Some(false));
        x.as_mut().map(|a| {a.pow(705737); a});
        assert_eq!(x.map(|a| a.is_one()), Some(true));
    }
}

