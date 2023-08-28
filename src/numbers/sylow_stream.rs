use crate::numbers::*;
use crate::util::*;

const STACK_ADDITION_LIMIT: u8 = 127;

pub mod flags {
    pub const NONE: u8 = 0x01;
    pub const NO_UPPER_HALF: u8 = 0x02;
    pub const LEQ: u8 = 0x04;
    pub const NO_PARABOLIC: u8 = 0x08;
}

pub struct SylowStreamBuilder<'a, C: SylowDecomposable + std::fmt::Debug> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,
    targets: Vec<Vec<u128>>,
    has_trivial: bool
}

pub struct SylowStream<'a, C: SylowDecomposable + std::fmt::Debug> {
    decomp: &'a SylowDecomp<'a, C>,
    mode: u8,
    targets: Vec<Vec<u128>>,
    stack: Vec<StackElem<'a, C>>,
}

#[derive(Debug)]
struct Seed<'a, C: SylowDecomposable> {
    i: usize,
    step: u128,
    rs: Vec<u128>,
    coords: SylowElem<'a, C>,
    block_upper: bool,
    contributed: bool,
}

#[derive(Debug)]
enum StackElem<'a, C: SylowDecomposable> {
    Res(SylowElem<'a, C>),
    Seed(Seed<'a, C>),
    Thunk {
        seed: Seed<'a, C>,
        start: u128,
        stop: u128,
    },
}

mod statuses {
    pub const EQ: u8 = 0x01;
    pub const ONE_AWAY: u8 = 0x02;
    pub const KEEP_GOING: u8 = 0x04;
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> SylowStreamBuilder<'a, C> {
    pub fn new(decomp: &'a SylowDecomp<'a, C>) -> SylowStreamBuilder<'a, C> {
        SylowStreamBuilder {
            decomp,
            mode: flags::NONE,
            targets: Vec::new(),
            has_trivial: false,
        }
    }

    pub fn add_flag(mut self, mode: u8) -> SylowStreamBuilder<'a, C> {
        self.mode |= mode;
        self
    }

    pub fn add_target(mut self, t: Vec<u128>) -> SylowStreamBuilder<'a, C> {
        if t.iter().all(|x| *x == 0) {
            self.has_trivial = true;
        } else {
            self.targets.push(t);
        }
        self
    }

    pub fn build(&self) -> SylowStream<'a, C> {
        let mut stack = Vec::new();

        if self.has_trivial || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0) {
            stack.push(StackElem::Res(self.decomp.one()));
        }

        for i in 0..self.decomp.factors().len() {
            let (p, d) = self.decomp.factors()[i];

            'a: for t in &self.targets {
                if self.mode & flags::LEQ == 0 {
                    for t in &t[0..i] {
                        if *t > 0 {
                            continue 'a;
                        }
                    }
                }

                let mut coords = self.decomp.one();
                let mut rs = vec![0; self.decomp.factors().len()];
                let mut step = intpow(p, d - 1, 0);

                if self.mode & flags::NO_PARABOLIC != 0 {
                    if t[i] == 0 {
                        continue;
                    }
                    if p == 2 && t[i] == 1 {
                        continue;
                    }
                    if p == 2 && t[i] > 1 {
                        coords.coords[i] += step;
                        rs[i] += 1;
                        step >>= 1;
                    }
                }

                stack.push(StackElem::Seed(Seed {
                    i,
                    coords,
                    step,
                    rs,
                    block_upper: self.mode & flags::NO_UPPER_HALF != 0,
                    contributed: false,
                }));
                break;
            }
        }
        SylowStream {
            decomp: self.decomp,
            mode: self.mode,
            targets: self.targets.clone(),
            stack,
        }
    }
}

impl<'a, C: SylowDecomposable + std::fmt::Debug> Iterator for SylowStream<'a, C> {
    type Item = SylowElem<'a, C>;

    fn next(&mut self) -> Option<SylowElem<'a, C>> {
        let Some(top) = self.stack.pop() else { return None; };
        // println!("top is {:?}", top);
        let (seed, pause) = match top {
            StackElem::Res(coords) => {
                return Some(coords);
            }
            StackElem::Seed(s) => (s, None),
            StackElem::Thunk { seed, start, stop } => (seed, Some((start, stop))),
        };

        let (p, _) = self.decomp.factors()[seed.i];

        let status = self.get_status(&seed.rs, seed.i);
        //println!("status: {status}");

        // First, create new seeds by incrementing
        // the current power.
        if status & statuses::KEEP_GOING != 0 {
            let (mut start, stop) = match pause {
                None => {
                    if status & statuses::ONE_AWAY != 0 && !self.has_flag(flags::LEQ) {
                        (1, p)
                    } else {
                        (0, p)
                    }
                }
                Some(x) => x,
            };

            if stop - start > STACK_ADDITION_LIMIT as u128 {
                self.stack.push(StackElem::Thunk {
                    seed: seed.clone(),
                    start,
                    stop: stop - STACK_ADDITION_LIMIT as u128,
                });
                start = stop - STACK_ADDITION_LIMIT as u128;
            }
            //println!("start: {start}");
            for j in start..stop {
                let mut coords = seed.coords.clone();
                coords.coords[seed.i] += j * seed.step;

                if seed.block_upper {
                    let fact = self.decomp.factors().factor(seed.i);
                    if coords.coords[seed.i] > fact / 2 {
                        break;
                    }
                }

                let mut rs = seed.rs.clone();
                rs[seed.i] += 1;

                let push = StackElem::Seed(Seed {
                    i: seed.i,
                    step: seed.step / p,
                    coords,
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
        if pause.is_none() && status & statuses::EQ != 0 && seed.contributed {
            let pushed = self.push_next_seeds(&seed, seed.i + 1);
            if self.has_flag(flags::LEQ) || !pushed {
                // println!("returning {:?}", seed.coords);
                return Some(seed.coords);
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
        let mut status = 0;
        for t in &self.targets {
            let skip = rs.iter().zip(t).take(i).any(|(r, t)| {
                self.has_flag(flags::LEQ) && r > t || !self.has_flag(flags::LEQ) && r != t
            });
            if skip {
                continue;
            }

            match t[i].overflowing_sub(rs[i]) {
                (0, false) => {
                    status |= statuses::EQ;
                }
                (1, false) => {
                    status |= statuses::ONE_AWAY | statuses::KEEP_GOING;
                }
                (_, false) => {
                    status |= statuses::KEEP_GOING;
                }
                (_, true) => {}
            }
        }
        status
    }

    fn push_next_seeds(&mut self, seed: &Seed<'a, C>, start: usize) -> bool {
        let mut pushed_any = false;
        // Note: In Rust, (a..a) is the empty iterator.
        for j in start..self.decomp.factors().len() {
            if !self.has_flag(flags::LEQ) {
                let status = self.get_status(&seed.rs, j);
                if status & statuses::KEEP_GOING == 0 {
                    continue;
                }
            }
            let (p, d) = self.decomp.factors()[j];
            let coords = seed.coords.clone();
            let s = Seed {
                i: j,
                coords,
                step: intpow(p, d - 1, 0),
                rs: seed.rs.clone(),
                block_upper: false,
                contributed: false,
            };
            // println!("pushing {:?}", s);
            self.stack.push(StackElem::Seed(s));
            pushed_any = true;
        }
        pushed_any
    }
}

impl<'a, C: SylowDecomposable> Clone for Seed<'a, C> {
    fn clone(&self) -> Seed<'a, C> {
        Seed {
            i: self.i,
            step: self.step,
            rs: self.rs.clone(),
            coords: self.coords.clone(),
            block_upper: self.block_upper,
            contributed: self.contributed,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::fp::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    pub fn test_make_stream() {
        let fp = Factorization::new(vec![(2, 1), (3, 1)]);
        let g = SylowDecomp::new(&FpStar::<7> {}, fp.clone());
        println!("g is {g:?}");
        let mut stream = SylowStreamBuilder::new(&g).add_target(vec![1, 0]).build();
        assert_eq!(
            stream.next().map(|s| {
                println!("s is {s:?}");
                s.to_product(&g).into()
            }),
            Some(6)
        );
        assert_eq!(stream.next(), None);
    }

    #[test]
    pub fn test_generates_small() {
        let fp = Factorization::new(vec![(2, 2), (3, 1), (5, 1)]);
        let g = SylowDecomp::new(&FpStar::<61> {}, fp.clone());
        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 0, 0])
            .build();
        let coords: Vec<SylowElem<FpStar<61>>> = stream.collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(!x.is_one(&g));
        x = x.pow(2, &g);
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
                x = x.multiply(&y, &g);
            }
            assert!(x.is_one(&g));
        }
        assert_eq!(count, 2);

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 1, 0])
            .build();
        let coords: Vec<SylowElem<FpStar<61>>> = stream.collect();
        assert_eq!(coords.len(), 2);
    }

    #[test]
    pub fn test_generates_big() {
        let fp = Factorization::new(vec![
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
        let g = SylowDecomp::new(&FpStar::<BIG_P> {}, fp.clone());

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 0, 0, 2, 0, 0, 0])
            .build();
        let coords: Vec<SylowElem<FpStar<BIG_P>>> = stream.collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        let mut stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 0, 0, 0, 0, 1, 0])
            .build();
        let mut x = stream.next();
        assert_eq!(x.as_ref().map(|a| a.is_one(&g)), Some(false));
        x = x.as_mut().map(|a| a.pow(705737, &g));
        assert_eq!(x.map(|a| a.is_one(&g)), Some(true));
    }

    #[test]
    pub fn test_generates_medium() {
        let fp = Factorization::new(vec![(2, 1), (3, 3), (5, 1)]);
        let g = SylowDecomp::new(&FpStar::<271> {}, fp.clone());

        let builder = SylowStreamBuilder::new(&g).add_target(vec![0, 2, 1]);
        let stream_all = builder.build();
        let coords: Vec<SylowElem<FpStar<271>>> = stream_all.collect();
        //println!("{coords:?}");
        assert_eq!(coords.len(), 24);
    }

    #[test]
    pub fn test_skips_upper_half() {
        let fp = Factorization::new(vec![(2, 1), (3, 3), (5, 1)]);
        let g = SylowDecomp::new(&FpStar::<271> {}, fp.clone());

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 2, 1])
            .add_flag(flags::NO_UPPER_HALF)
            .build();
        let coords: Vec<SylowElem<FpStar<271>>> = stream.collect();
        //println!("{coords:?}");
        assert_eq!(coords.len(), 12);
    }

    #[test]
    pub fn test_multiple_targets() {
        let fp = Factorization::new(vec![(2, 1), (3, 3), (5, 1)]);
        let g = SylowDecomp::new(&FpStar::<271> {}, fp.clone());

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![1, 0, 0])
            .add_target(vec![0, 1, 0])
            .build();
        let coords: Vec<SylowElem<FpStar<271>>> = stream.collect();
        assert_eq!(coords.len(), 3);

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
        let coords: Vec<SylowElem<FpStar<271>>> = stream.collect();
        for x in &coords {
            println!("{x:?}");
        }
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2() {
        let fp = Factorization::new(vec![(2, 1), (7, 2), (13, 2), (29, 2)]);
        let g = SylowDecomp::new(&FpStar::<13928643> {}, fp.clone());

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 1, 1, 0])
            .add_flag(flags::LEQ)
            .build();
        let coords: Vec<SylowElem<FpStar<13928643>>> = stream.collect();
        for x in &coords {
            println!("{x:?}");
        }
        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic() {
        let fact = Factorization::new(vec![(2, 2), (3, 1), (5, 1)]);
        let g = SylowDecomp::new(&FpStar::<61> {}, fact.clone());

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![2, 1, 0])
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .build();
        for mut x in stream {
            assert!(!x.is_one(&g));
            x = x.square(&g);
            assert!(!x.is_one(&g));
        }
    }
}
