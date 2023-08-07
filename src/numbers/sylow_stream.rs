use crate::numbers::sylow::*;
use crate::util::*;
use crate::numbers::factorization::*;

pub mod flags {
    pub const NONE: u8 = 0x01;
    pub const NO_UPPER_HALF: u8 = 0x02;
    pub const LEQ: u8 = 0x04;
}

#[derive(Debug)]
struct Seed<'a, C: SylowDecomposable> {
    i: usize,
    step: u128,
    rs: Vec<u128>,
    res: SylowElem<'a, C>,
    block_upper: bool,
    contributed: bool
}

#[derive(Debug)]
enum StackElem<'a, C: SylowDecomposable> {
    Res(SylowElem<'a, C>),
    Seed(Seed<'a, C>)
}

pub struct SylowStreamBuilder<'a, C: SylowDecomposable + std::fmt::Debug> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,
    targets: Vec<Vec<u128>>
}

pub struct SylowStream<'a, C: SylowDecomposable + std::fmt::Debug> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,
    targets: Vec<Vec<u128>>,
    stack: Vec<StackElem<'a, C>>
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> SylowStreamBuilder<'a, C> {
    pub fn new(decomp: &'a SylowDecomp<'a, C>) -> SylowStreamBuilder<'a, C> {
        SylowStreamBuilder {
            decomp,
            mode: flags::NONE,
            targets: Vec::new()
        }
    }

    pub fn add_flag(mut self, mode: u8) -> SylowStreamBuilder<'a, C> {
        self.mode |= mode;
        self
    }

    pub fn add_target(mut self, t: Vec<u128>) -> SylowStreamBuilder<'a, C> {
        self.targets.push(t);
        self
    }

    pub fn build(&self) -> SylowStream<'a, C> {
        let mut stack = Vec::new();
        if self.mode & flags::LEQ != 0 {
            stack.push(StackElem::Res(self.decomp.one()));
        }
        for i in 0..self.decomp.factors().len() {
            let (p,d) = self.decomp.factors()[i];
            'a: for t in &self.targets {
                if self.mode & flags::LEQ == 0 {
                    for t in &t[0..i] {
                        if *t > 0 { continue 'a; }
                    }
                }
                if t[i] == 0 { continue; }
                stack.push(StackElem::Seed(Seed {
                    i,
                    res: self.decomp.one(),
                    step: intpow(p, d - 1, 0),
                    rs: vec![0 ; self.decomp.factors().len()],
                    block_upper: self.mode & flags::NO_UPPER_HALF != 0,
                    contributed: false
                }));
                break;
            }
        }
        SylowStream {
            decomp: self.decomp,
            mode: self.mode,
            targets: self.targets.clone(),
            stack
        }
    }
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> Iterator for SylowStream<'a, C> {
    type Item = SylowElem<'a, C>;

    fn next(&mut self) -> Option<SylowElem<'a, C>> {
        let Some(top) = self.stack.pop() else { return None; };
        // println!("top is {:?}", top);
        let seed = match top {
            StackElem::Res(res) => { 
                return Some(res); 
            }
            StackElem::Seed(s) => { s } 
        };

        let (p,_) = self.decomp.factors()[seed.i];

        let status = self.get_status(&seed.rs, seed.i);
        //println!("status: {status}");

        // First, create new seeds by incrementing
        // the current power.
        if status & statuses::KEEP_GOING != 0 {
            let start = if status & statuses::ONE_AWAY != 0 
                && !self.has_flag(flags::LEQ) 
                { 1 } else { 0 };
            //println!("start: {start}");
            for j in start..p {
                let mut res = seed.res.clone();
                res.coords[seed.i] += j * seed.step;

                if seed.block_upper {
                    let fact = self.decomp.factors().factor(seed.i);
                    if res.coords[seed.i] > fact / 2 { break; }
                }

                let mut rs = seed.rs.clone();
                rs[seed.i] += 1;

                let push = StackElem::Seed(Seed {
                    i: seed.i,
                    step: seed.step / p,
                    res,
                    rs,
                    block_upper: seed.block_upper,
                    contributed: j != 0,
                });
                // println!("pushing {:?}", push);
                self.stack.push(push);
            }
        }

        // Next, create new seeds by moving to the next prime power,
        // but only if we are *done* with this prime power.
        if status & statuses::EQ != 0 && seed.contributed {
            let pushed = self.push_next_seeds(&seed, seed.i + 1);
            if self.has_flag(flags::LEQ) || !pushed {
                // println!("returning {:?}", seed.res);
                return Some(seed.res);
            }
        }
        self.next()
    }
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> SylowStream<'a, C> {
    fn has_flag(&self, flag: u8) -> bool {
        self.mode & flag != 0
    }

    fn get_status(&self, rs: &[u128], i: usize) -> u8 {
        let mut res = 0;
        for t in &self.targets {
            let skip = rs
                .iter()
                .zip(t)
                .take(i)
                .any(|(r, t)| {
                    self.has_flag(flags::LEQ) && r > t
                    || !self.has_flag(flags::LEQ) && r != t
                });
            if skip { continue; }

            match t[i].overflowing_sub(rs[i]) {
                (0,false) => { res |= statuses::EQ; }
                (1,false) => { res |= statuses::ONE_AWAY | statuses::KEEP_GOING; }
                (_,false) => { res |= statuses::KEEP_GOING; }
                (_,true) => {}
            }
        }
        res
    }

    fn push_next_seeds(&mut self, seed: &Seed<'a, C>, start: usize) -> bool {
        let mut pushed_any = false;
        // Note: In Rust, (a..a) is the empty iterator.
        for j in start..self.decomp.factors().len() {
            if !self.has_flag(flags::LEQ) {
                let status = self.get_status(&seed.rs, j);
                if status & statuses::KEEP_GOING == 0 { continue; }
            }
            let (p,d) = self.decomp.factors()[j];
            let res = seed.res.clone();
            let s = Seed {
                i: j,
                res,
                step: intpow(p, d - 1, 0),
                rs: seed.rs.clone(),
                block_upper: false,
                contributed: false
            };
            // println!("pushing {:?}", s);
            self.stack.push(StackElem::Seed(s));
            pushed_any = true;
        }
        pushed_any
    }
}

mod statuses {
    pub const EQ: u8 = 0x01;
    pub const ONE_AWAY: u8 = 0x02;
    pub const KEEP_GOING: u8 = 0x04;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::fp::*;

    #[test]
    pub fn test_make_stream() {
        let fp = Factorization::new(vec![(2, 1), (3, 1)]);
        let g = SylowDecomp::new(&fp);
        let mut stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 0])
            .build();
        assert_eq!(stream.next().map(|s| s.to_product(&g).value()), Some(6));
        assert_eq!(stream.next(), None);
    }

    #[test]
    pub fn test_generates_small() {
        let fp = Factorization::new(vec![(2, 2), (3, 1), (5, 1)]);
        let g = SylowDecomp::new(&fp);
        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 0, 0]) 
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        assert_eq!(res.len(), 1);
        let mut x = res[0].clone();
        assert!(!x.is_one(&g));
        x.pow(2, &g);
        assert!(x.is_one(&g));

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![2, 0, 0]) 
            .build();
        let mut count = 0;
        for mut x in stream {
            //println!("streamed {:?}", x);
            count += 1;
            for _ in 1..3 {
                let y = x.clone();
                assert!(!x.is_one(&g));
                x.multiply(&y, &g);
            }
            assert!(x.is_one(&g));
        }
        assert_eq!(count, 2);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 1, 0])
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        assert_eq!(res.len(), 2);
    }

    #[test]
    pub fn test_generates_big() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let g = SylowDecomp::new(&fp);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 0, 0, 2, 0, 0, 0])
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        assert_eq!(res.len(), 29 * 29 - 29);

        let mut stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 0, 0, 0, 0, 1, 0])
            .build();
        let mut x = stream.next();
        assert_eq!(x.as_ref().map(|a| a.is_one(&g)), Some(false));
        x.as_mut().map(|a| {a.pow(705737, &g); a});
        assert_eq!(x.map(|a| a.is_one(&g)), Some(true));
    }

    #[test]
    pub fn test_generates_medium() {
        let fp = Factorization::new(vec![(2,1), (3,3), (5,1)]);
        let g = SylowDecomp::new(&fp);

        let builder = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 2, 1]);
        let stream_all = builder.build();
        let res: Vec<SylowElem<FpStar>> = stream_all.collect();
        //println!("{res:?}");
        assert_eq!(res.len(), 24);
    }

    #[test]
    pub fn test_skips_upper_half() {
        let fp = Factorization::new(vec![(2,1), (3,3), (5,1)]);
        let g = SylowDecomp::new(&fp);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 2, 1])
            .add_flag(flags::NO_UPPER_HALF)
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        //println!("{res:?}");
        assert_eq!(res.len(), 12);
    }

    #[test]
    pub fn test_multiple_targets() {
        let fp = Factorization::new(vec![(2, 1), (3, 3), (5, 1)]);
        let g = SylowDecomp::new(&fp);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 0, 0])
            .add_target(vec![0, 1, 0])
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        assert_eq!(res.len(), 3);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 1, 0])
            .add_target(vec![0, 2, 0])
            .add_target(vec![0, 0, 1])
            .add_flag(flags::LEQ)
            .build();
        /*
        for x in stream {
            println!("{x:?}");
        }
        assert!(false);
        */
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        for x in &res {
            println!("{x:?}");
        }
        assert_eq!(res.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2() {
        let fp = Factorization::new(vec![(2, 1), (7, 2), (13, 2), (29, 2)]);
        let g = SylowDecomp::new(&fp);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 1, 1, 0])
            .add_flag(flags::LEQ)
            .build();
        let res: Vec<SylowElem<FpStar>> = stream.collect();
        for x in &res {
            println!("{x:?}");
        }
        assert_eq!(res.len(), 91);
    }
}


