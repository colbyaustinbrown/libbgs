use crossbeam::deque::{Injector, Stealer, Worker};
use crossbeam::queue::SegQueue;
use rayon::iter::plumbing;
use rayon::iter::plumbing::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use std::cell::RefCell;

use crate::numbers::*;
use crate::util::*;

const STACK_ADDITION_LIMIT: u8 = 127;

/// Bitwise flags for configuring a SylowStreamBuilder.
/// May be combined with the bitwise integer OR operator, `|`.
pub mod flags {
    /// Stream will behave with all default options.
    /// Identity when used with the `|` operator.
    pub const NONE: u8 = 0x01;

    /// Stream will yield half as many elements as the default.
    /// Precisely, the stream will yield either $\chi$ or $\chi^{-1}$, but not both, for every
    /// $\chi \in F_{p^2}$.
    /// Even more precisely, for every $x$ the builder yields, the first non-zero coordinate is
    /// guarantied to be less than half the maximum of the its corresponding prime power.
    pub const NO_UPPER_HALF: u8 = 0x02;

    /// Stream will yield any element less than or equal to a target.
    /// Guaranteed to only  return 1 value of `\chi` for each threat; guaranteed 1 `for` 1.
    pub const LEQ: u8 = 0x04;

    ///  Stream will yield elements in a parabolic order.
    pub const NO_PARABOLIC: u8 = 0x08;

    /// Stream is guaranteed to include the identity.
    pub const INCLUDE_ONE: u8 = 0x10;
}

/// A builder for a stream yielding elements of particular orders, as their Sylow decompositions.
pub struct SylowStreamBuilder<'a, S, const L: usize, C: SylowDecomposable<S, L> + std::fmt::Debug> {
    decomp: &'a SylowDecomp<S, L, C>,
    mode: u8,
    targets: Vec<[u128; L]>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements in parallel on multiple threads.
pub struct SylowParStream<
    'a,
    S: Send + Sync,
    const L: usize,
    C: SylowDecomposable<S, L> + std::fmt::Debug,
> {
    decomp: &'a SylowDecomp<S, L, C>,
    mode: u8,
    targets: Vec<[u128; L]>,
    global_stack: Injector<Seed<S, L, C>>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements sequentially on a single thread.
pub struct SylowSeqStream<'a, S, const L: usize, C: SylowDecomposable<S, L>> {
    decomp: &'a SylowDecomp<S, L, C>,
    mode: u8,
    targets: Vec<[u128; L]>,
    stack: Vec<Seed<S, L, C>>,
    buffer: Vec<SylowElem<S, L, C>>,
}

struct SylowStreamWorker<'a, S, const L: usize, C, Con>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L>,
    Con: Consumer<SylowElem<S, L, C>>,
{
    stream: &'a SylowParStream<'a, S, L, C>,
    local_stack: Worker<Seed<S, L, C>>,
    stealers: &'a [Stealer<Seed<S, L, C>>],
    folder: RefCell<Option<Con::Folder>>,
}

#[derive(Debug)]
struct Seed<S, const L: usize, C: SylowDecomposable<S, L>> {
    i: usize,
    step: u128,
    rs: [u128; L],
    part: SylowElem<S, L, C>,
    block_upper: bool,
    start: u128,
}

struct Status(u8);

mod statuses {
    pub const EQ: u8 = 0x01;
    pub const ONE_AWAY: u8 = 0x02;
    pub const KEEP_GOING: u8 = 0x04;
}

trait SylowStream<'a, S, const L: usize, C>
where
    C: SylowDecomposable<S, L>,
{
    fn decomp(&self) -> &'a SylowDecomp<S, L, C>;
    fn mode(&self) -> u8;
    fn targets(&self) -> &[[u128; L]];
    fn push(&mut self, e: Seed<S, L, C>);
    fn consume(&mut self, e: SylowElem<S, L, C>);

    fn has_flag(&self, flag: u8) -> bool {
        self.mode() & flag != 0
    }

    fn get_start(&self, status: &Status) -> u128 {
        if status.0 & statuses::ONE_AWAY != 0 && !self.has_flag(flags::LEQ) {
            1
        } else {
            0
        }
    }

    fn propogate(&mut self, seed: Seed<S, L, C>) {
        let (p, _) = <C as Factored<S, L>>::FACTORS[seed.i];

        let status = self.get_status(&seed.rs, seed.i);
        //println!("status: {status}");

        // First, create new seeds by incrementing
        // the current power.
        if status.has(statuses::KEEP_GOING) {
            let mut stop = p;

            if stop - seed.start > STACK_ADDITION_LIMIT as u128 {
                self.push(Seed {
                    start: seed.start + STACK_ADDITION_LIMIT as u128,
                    ..seed
                });
                stop = seed.start + STACK_ADDITION_LIMIT as u128;
            }
            let lim = <C as Factored<S, L>>::FACTORS.factor(seed.i) / 2;
            for j in seed.start..stop {
                let tmp = seed.part.coords[seed.i] + j * seed.step;

                if seed.block_upper && tmp > lim {
                    break;
                }

                let mut part = seed.part;
                part.coords[seed.i] = tmp;

                let mut rs = seed.rs;
                rs[seed.i] += 1;

                let status = self.get_status(&rs, seed.i);
                let next = Seed {
                    i: seed.i,
                    step: seed.step / p,
                    part,
                    rs,
                    block_upper: seed.block_upper,
                    start: self.get_start(&status),
                };

                // Next, create new seeds by moving to the next prime power,
                // but only if we are *done* with this prime power.
                if status.has(statuses::EQ) && j != 0 {
                    let mut pushed_any = false;
                    // Note: In Rust, (a..a) is the empty iterator.
                    for k in (next.i + 1)..L {
                        let status = self.get_status(&next.rs, k);
                        if !(self.has_flag(flags::LEQ) || status.has(statuses::KEEP_GOING)) {
                            continue;
                        }
                        let (p, d) = <C as Factored<S, L>>::FACTORS[k];
                        let s = Seed {
                            i: k,
                            part: next.part,
                            step: intpow(p, d - 1, 0),
                            rs: next.rs,
                            block_upper: false,
                            start: self.get_start(&status),
                        };
                        self.push(s);
                        pushed_any = true;
                    }
                    if self.has_flag(flags::LEQ) || !pushed_any {
                        self.consume(next.part);
                    }
                }
                if status.has(statuses::KEEP_GOING) {
                    self.push(next);
                }
            }
        }
    }

    fn get_status(&self, rs: &[u128], i: usize) -> Status {
        let mut status = 0;
        for t in self.targets() {
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
        Status(status)
    }
}

impl<'a, S, const L: usize, C: SylowDecomposable<S, L> + std::fmt::Debug>
    SylowStreamBuilder<'a, S, L, C>
{
    fn get_starting_stack(&self) -> Vec<Seed<S, L, C>> {
        let mut res = Vec::new();

        for i in 0..L {
            let (p, d) = <C as Factored<S, L>>::FACTORS[i];

            'a: for t in &self.targets {
                if self.mode & flags::LEQ == 0 {
                    for t in &t[0..i] {
                        if *t > 0 {
                            continue 'a;
                        }
                    }
                }

                let mut part = SylowElem::one();
                let mut rs = [0; L];
                let mut step = intpow(p, d - 1, 0);

                if self.mode & flags::NO_PARABOLIC != 0 {
                    if t[i] == 0 {
                        continue;
                    }
                    if p == 2 && t[i] == 1 {
                        continue;
                    }
                    if p == 2 && t[i] > 1 {
                        part.coords[i] += step;
                        rs[i] += 1;
                        step >>= 1;
                    }
                }

                let start =
                    if t[i].overflowing_sub(rs[i]) == (1, false) && self.mode & flags::LEQ == 0 {
                        1
                    } else {
                        0
                    };
                res.push(Seed {
                    i,
                    part,
                    step,
                    rs,
                    block_upper: self.mode & flags::NO_UPPER_HALF != 0,
                    start,
                });
                break;
            }
        }
        res
    }

    /// Returns a new `SylowStreamBuilder`, yielding elements of the type associated with `decomp`.
    pub fn new(decomp: &'a SylowDecomp<S, L, C>) -> SylowStreamBuilder<'a, S, L, C> {
        SylowStreamBuilder {
            decomp,
            mode: flags::NONE,
            targets: Vec::new(),
        }
    }

    /// Adds a flag to the `SylowStreamBuilder`, modifying its yields.
    pub fn add_flag(mut self, mode: u8) -> SylowStreamBuilder<'a, S, L, C> {
        self.mode |= mode;
        self
    }

    /// Adds a target order to this `SylowStreamBuilder`.
    /// The `SylowStream` built from this builder will only yield elements of the orders of
    /// `target`s, or elements of order dividing `target` if `target
    pub fn add_target(mut self, t: [u128; L]) -> SylowStreamBuilder<'a, S, L, C> {
        if t.iter().all(|x| *x == 0) {
            self.mode |= flags::INCLUDE_ONE;
        } else {
            self.targets.push(t);
        }
        self
    }
}

impl<'a, S, const L: usize, C> ParallelIterator for SylowParStream<'a, S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L> + Send + Sync,
{
    type Item = SylowElem<S, L, C>;

    fn drive_unindexed<Con>(self, consumer: Con) -> Con::Result
    where
        Con: plumbing::UnindexedConsumer<Self::Item>,
    {
        let consumer_spawner = consumer.split_off_left();

        let num_threads = std::thread::available_parallelism().unwrap().into();
        let mut consumers = Vec::new();
        let mut workers: Vec<Worker<Seed<S, L, C>>> = Vec::new();
        let mut stealers = Vec::new();

        for _ in 0..num_threads {
            let con = consumer_spawner.split_off_left();
            let wrk = Worker::new_fifo();
            let stl = wrk.stealer();

            consumers.push(con);
            workers.push(wrk);
            stealers.push(stl);
        }

        let results = SegQueue::new();
        let mut folder = consumer.into_folder();
        if self.has_flag(flags::INCLUDE_ONE)
            || (self.has_flag(flags::LEQ) && !self.has_flag(flags::NO_PARABOLIC))
        {
            folder = folder.consume(SylowElem::one());
            results.push(folder.complete());
        }
        rayon::scope(|s| {
            workers
                .into_iter()
                .zip(consumers.into_iter())
                .for_each(|(wrk, consumer)| {
                    s.spawn(|_| {
                        let folder = consumer.into_folder();
                        let worker: SylowStreamWorker<S, L, C, Con> = SylowStreamWorker {
                            stream: &self,
                            local_stack: wrk,
                            stealers: &stealers,
                            folder: RefCell::new(Some(folder)),
                        };
                        results.push(worker.work());
                    });
                });
        });

        let mut res = results.pop().unwrap();
        for r in results {
            res = consumer_spawner.to_reducer().reduce(res, r);
        }
        res
    }
}

impl<
        'a,
        S: Send + Sync,
        const L: usize,
        C: SylowDecomposable<S, L>,
        Con: Consumer<SylowElem<S, L, C>>,
    > SylowStreamWorker<'a, S, L, C, Con>
{
    fn work(mut self) -> Con::Result {
        let mut retry = 0;
        loop {
            while let Some(task) = self.local_stack.pop().or_else(|| {
                std::iter::repeat_with(|| {
                    self.stream
                        .global_stack
                        .steal_batch_and_pop(&self.local_stack)
                        .or_else(|| self.stealers.iter().map(|s| s.steal()).collect())
                })
                .find(|s| !s.is_retry())
                .and_then(|s| s.success())
            }) {
                self.propogate(task);
            }
            if retry < 3 {
                std::thread::sleep(std::time::Duration::from_millis(5));
                retry += 1;
                continue;
            }
            break;
        }
        self.folder.into_inner().unwrap().complete()
    }
}

impl<'a, S, const L: usize, C: SylowDecomposable<S, L>> Iterator for SylowSeqStream<'a, S, L, C> {
    type Item = SylowElem<S, L, C>;

    fn next(&mut self) -> Option<SylowElem<S, L, C>> {
        if let Some(res) = self.buffer.pop() {
            Some(res)
        } else if let Some(top) = self.stack.pop() {
            self.propogate(top);
            self.next()
        } else {
            None
        }
    }
}

impl<'a, S, const L: usize, C> IntoParallelIterator for SylowStreamBuilder<'a, S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L> + Send + Sync,
{
    type Item = SylowElem<S, L, C>;
    type Iter = SylowParStream<'a, S, L, C>;

    /// Returns a SylowStream yielding the elements requested via the constructor, `add_flag`, and
    /// `flag_target` invocations.
    fn into_par_iter(self) -> SylowParStream<'a, S, L, C> {
        let global_stack = Injector::new();

        self.get_starting_stack()
            .into_iter()
            .for_each(|x| global_stack.push(x));

        SylowParStream {
            decomp: self.decomp,
            mode: self.mode,
            targets: self.targets.clone(),
            global_stack,
        }
    }
}

impl<'a, S, const L: usize, C> IntoIterator for SylowStreamBuilder<'a, S, L, C>
where
    C: SylowDecomposable<S, L>,
{
    type Item = SylowElem<S, L, C>;
    type IntoIter = SylowSeqStream<'a, S, L, C>;

    fn into_iter(self) -> SylowSeqStream<'a, S, L, C> {
        let mut stream = SylowSeqStream {
            decomp: self.decomp,
            mode: self.mode,
            targets: self.targets.clone(),
            stack: self.get_starting_stack(),
            buffer: Vec::new(),
        };
        if (stream.mode & flags::INCLUDE_ONE != 0)
            || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0)
        {
            stream.buffer.push(SylowElem::one());
        }
        stream
    }
}

impl Status {
    fn has(&self, flag: u8) -> bool {
        self.0 & flag != 0
    }
}

impl<'a, S, const L: usize, C> SylowStream<'a, S, L, C> for SylowParStream<'a, S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L>,
{
    fn push(&mut self, e: Seed<S, L, C>) {
        self.global_stack.push(e);
    }

    fn consume(&mut self, _: SylowElem<S, L, C>) {
        unimplemented!();
    }

    fn decomp(&self) -> &'a SylowDecomp<S, L, C> {
        self.decomp
    }

    fn mode(&self) -> u8 {
        self.mode
    }

    fn targets(&self) -> &[[u128; L]] {
        &self.targets
    }
}

impl<'a, S, const L: usize, C, Con> SylowStream<'a, S, L, C> for SylowStreamWorker<'a, S, L, C, Con>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L>,
    Con: Consumer<SylowElem<S, L, C>>,
{
    fn push(&mut self, e: Seed<S, L, C>) {
        self.local_stack.push(e);
    }

    fn consume(&mut self, e: SylowElem<S, L, C>) {
        let mut f = self.folder.take().unwrap();
        f = f.consume(e);
        self.folder.replace(Some(f));
    }

    fn decomp(&self) -> &'a SylowDecomp<S, L, C> {
        self.stream.decomp
    }

    fn mode(&self) -> u8 {
        self.stream.mode
    }

    fn targets(&self) -> &[[u128; L]] {
        &self.stream.targets
    }
}

impl<'a, S, const L: usize, C> SylowStream<'a, S, L, C> for SylowSeqStream<'a, S, L, C>
where
    C: SylowDecomposable<S, L>,
{
    fn push(&mut self, e: Seed<S, L, C>) {
        self.stack.push(e);
    }

    fn consume(&mut self, e: SylowElem<S, L, C>) {
        self.buffer.push(e);
    }

    fn decomp(&self) -> &'a SylowDecomp<S, L, C> {
        self.decomp
    }

    fn mode(&self) -> u8 {
        self.mode
    }

    fn targets(&self) -> &[[u128; L]] {
        &self.targets
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Clone for Seed<S, L, C> {
    fn clone(&self) -> Seed<S, L, C> {
        Seed {
            i: self.i,
            step: self.step,
            rs: self.rs,
            part: self.part,
            block_upper: self.block_upper,
            start: self.start,
        }
    }
}
impl<S, const L: usize, C: SylowDecomposable<S, L>> Copy for Seed<S, L, C> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::fp::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    impl Factored<Phantom, 2> for FpNum<7> {
        const FACTORS: Factorization<2> = Factorization::new([(2, 1), (3, 1)]);
    }

    impl Factored<Phantom, 3> for FpNum<61> {
        const FACTORS: Factorization<3> = Factorization::new([(2, 2), (3, 1), (5, 1)]);
    }

    impl Factored<Phantom, 7> for FpNum<BIG_P> {
        const FACTORS: Factorization<7> = Factorization::new([
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
    }

    impl Factored<Phantom, 3> for FpNum<271> {
        const FACTORS: Factorization<3> = Factorization::new([(2, 1), (3, 3), (5, 1)]);
    }

    impl Factored<Phantom, 4> for FpNum<13928643> {
        const FACTORS: Factorization<4> = Factorization::new([(2, 1), (7, 2), (13, 2), (29, 2)]);
    }

    #[test]
    pub fn test_make_stream_seq() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new(&g)
            .add_target([1, 0])
            .into_iter()
            .filter_map(|s| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_seq() {
        let g = SylowDecomp::new();
        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 0, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<61>>> = stream.collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(!x.is_one());
        x = x.pow(2);
        assert!(x.is_one());

        let mut count = 0;
        SylowStreamBuilder::new(&g)
            .add_target([2, 0, 0])
            .into_iter()
            .for_each(|mut x| {
                count += 1;
                for _ in 1..3 {
                    let y = x.clone();
                    assert!(!x.is_one());
                    x = x.multiply(&y);
                }
                assert!(x.is_one());
            });
        assert_eq!(count, 2);

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 1, 0])
            .into_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_generates_big_seq() {
        let g = SylowDecomp::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 0, 0, 2, 0, 0, 0])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::new(&g)
            .add_target([0, 0, 0, 0, 0, 1, 0])
            .into_iter()
            .take(2)
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.pow(705737);
                assert!(x.is_one());
            });
    }

    #[test]
    pub fn test_generates_medium_seq() {
        let g: SylowDecomp<Phantom, 3, FpNum<271>> = SylowDecomp::new();
        let builder = SylowStreamBuilder::new(&g).add_target([0, 2, 1]);
        let stream_all = builder.into_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_seq() {
        let g = SylowDecomp::<Phantom, 3, FpNum<271>>::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 2, 1])
            .add_flag(flags::NO_UPPER_HALF)
            .into_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_seq() {
        let g = SylowDecomp::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 0, 0])
            .add_target([0, 1, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 1, 0])
            .add_target([0, 2, 0])
            .add_target([0, 0, 1])
            .add_flag(flags::LEQ)
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_seq() {
        let g = SylowDecomp::new();

        let coords: Vec<SylowElem<Phantom, 4, FpNum<13928643>>> = SylowStreamBuilder::new(&g)
            .add_target([0, 1, 1, 0])
            .add_flag(flags::LEQ)
            .into_iter()
            .collect();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_seq() {
        let g = SylowDecomp::<Phantom, 3, FpNum<61>>::new();

        SylowStreamBuilder::new(&g)
            .add_target([2, 1, 0])
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .into_iter()
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.square();
                assert!(!x.is_one());
            });
    }
    #[test]
    pub fn test_make_stream_par() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new(&g)
            .add_target([1, 0])
            .into_par_iter()
            .filter_map(|s| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_par() {
        let g = SylowDecomp::new();
        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 0, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<61>>> = stream.collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(!x.is_one());
        x = x.pow(2);
        assert!(x.is_one());

        let count = AtomicUsize::new(0);
        SylowStreamBuilder::new(&g)
            .add_target([2, 0, 0])
            .into_par_iter()
            .for_each(|mut x| {
                count.fetch_add(1, Ordering::Relaxed);
                for _ in 1..3 {
                    let y = x.clone();
                    assert!(!x.is_one());
                    x = x.multiply(&y);
                }
                assert!(x.is_one());
            });
        assert_eq!(count.into_inner(), 2);

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 1, 0])
            .into_par_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_generates_big_par() {
        let g = SylowDecomp::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 0, 0, 2, 0, 0, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::new(&g)
            .add_target([0, 0, 0, 0, 0, 1, 0])
            .into_par_iter()
            .take_any(2)
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.pow(705737);
                assert!(x.is_one());
            });
    }

    #[test]
    pub fn test_generates_medium_par() {
        let g: SylowDecomp<Phantom, 3, FpNum<271>> = SylowDecomp::new();
        let builder = SylowStreamBuilder::new(&g).add_target([0, 2, 1]);
        let stream_all = builder.into_par_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_par() {
        let g = SylowDecomp::<Phantom, 3, FpNum<271>>::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([0, 2, 1])
            .add_flag(flags::NO_UPPER_HALF)
            .into_par_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_par() {
        let g = SylowDecomp::new();

        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 0, 0])
            .add_target([0, 1, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::new(&g)
            .add_target([1, 1, 0])
            .add_target([0, 2, 0])
            .add_target([0, 0, 1])
            .add_flag(flags::LEQ)
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_par() {
        let g = SylowDecomp::new();

        let coords: Vec<SylowElem<Phantom, 4, FpNum<13928643>>> = SylowStreamBuilder::new(&g)
            .add_target([0, 1, 1, 0])
            .add_flag(flags::LEQ)
            .into_par_iter()
            .collect();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_par() {
        let g = SylowDecomp::<Phantom, 3, FpNum<61>>::new();

        SylowStreamBuilder::new(&g)
            .add_target([2, 1, 0])
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .into_par_iter()
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.square();
                assert!(!x.is_one());
            });
    }
}
