use crate::sylow::*;
use crate::util::*;
use crate::factorization::*;

pub struct SylowFactory<'a, C: SylowDecomposable> {
    pub decomp: &'a SylowDecomp<'a, C>,
    pub i: usize,

    stack: Vec<(u128, u128, usize)>
}

impl<'a, C: SylowDecomposable> SylowFactory<'a, C> {
    pub fn new(decomp: &'a SylowDecomp<C>, i: usize, r: usize) -> SylowFactory<'a, C> {
        let (p,d) = decomp.factors()[i];
        let step = intpow(p, d - 1, 0);
        SylowFactory {
            decomp: decomp,
            i,
            stack: vec![(0,step,r)]
        }
    }
}

impl<'a, C: SylowDecomposable> Iterator for SylowFactory<'a, C> {
    type Item = SylowElem<'a, C>;

    fn next(&mut self) -> Option<SylowElem<'a, C>> {
        loop {
            let l = self.stack.len();
            if l == 0 { return None; }

            let p = self.decomp.factors()[self.i].0;
            let (x, step, r) = self.stack.remove(l - 1);
            if r == 0 {
                let mut coords = vec![0 ; self.decomp.factors().len()];
                coords[self.i] = x;
                return Some(SylowElem::new(coords));
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
    use crate::fp::*;

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

