use rayon::iter::plumbing::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use std::cell::RefCell;
use std::marker::PhantomData;
use std::sync::Arc;

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
pub struct SylowStreamBuilder<S, const L: usize, C: SylowDecomposable<S, L> + std::fmt::Debug> {
    mode: u8,
    tree: FactorNode<L>,
    _phantom: PhantomData<(S, C)>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements in parallel on multiple threads.
pub struct SylowParStream<
    S: Send + Sync,
    const L: usize,
    C: SylowDecomposable<S, L> + std::fmt::Debug,
> {
    stack: Vec<Seed<S, L, C>>,
    splits: usize,
    emit_one: bool,
    tree: Arc<FactorNode<L>>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements sequentially on a single thread.
pub struct SylowSeqStream<S, const L: usize, C: SylowDecomposable<S, L>> {
    stack: Vec<Seed<S, L, C>>,
    buffer: Vec<SylowElem<S, L, C>>,
    tree: FactorNode<L>,
}

#[derive(Debug)]
struct Seed<S, const L: usize, C: SylowDecomposable<S, L>> {
    part: SylowElem<S, L, C>,
    start: u128,
    node: *const FactorNode<L>,
}

#[derive(Clone, Debug)]
struct FactorNode<const L: usize> {
    i: usize,
    rs: [usize; L],
    consume: bool,
    step: u128,
    lim: u128,
    next: [Option<Box<FactorNode<L>>>; L],
}

trait SylowStream<'a, S, const L: usize, C>
where
    Self: Sized,
    C: SylowDecomposable<S, L>,
{
    fn push(&mut self, e: Seed<S, L, C>);
    fn stack(&mut self) -> &mut Vec<Seed<S, L, C>>;
    fn tree(&self) -> &FactorNode<L>;

    fn init_stack(&mut self, mode: u8) {
        for i in 0..L {
            let Some(n) = self.tree().next[i].as_ref() else {
                continue;
            };

            let seed = Seed {
                part: SylowElem::one(),
                start: 0,
                node: &**n,
            };

            let (p, _) = <C as Factor<S, L>>::FACTORS[i];
            if mode & flags::NO_PARABOLIC != 0 && p == 2 {
                self.propogate(seed, |_, _| {});
            } else {
                self.push(seed);
            }
        }
    }

    fn propogate<F>(&mut self, seed: Seed<S, L, C>, mut consume: F)
    where
        Self: Sized,
        F: FnMut(&mut Self, SylowElem<S, L, C>),
    {
        let node = unsafe { &*seed.node };
        let (p, _) = <C as Factor<S, L>>::FACTORS[node.i];

        // First, create new seeds by incrementing
        // the current power.
        let mut stop = p;

        if stop - seed.start > STACK_ADDITION_LIMIT as u128 {
            self.push(Seed {
                start: seed.start + STACK_ADDITION_LIMIT as u128,
                ..seed
            });
            stop = seed.start + STACK_ADDITION_LIMIT as u128;
        }

        for j in seed.start..stop {
            let tmp = seed.part.coords[node.i] + j * node.step;
            if tmp > node.lim {
                break;
            }
            let mut part = seed.part;
            part.coords[node.i] = tmp;

            if let Some(n) = node.next[node.i].as_ref() {
                self.push(Seed {
                    part,
                    start: 0,
                    node: &**n,
                });
            }

            // Next, create new seeds by moving to the next prime power,
            // but only if we are *done* with this prime power.
            if j == 0 {
                continue;
            }
            if node.consume {
                consume(self, part);
            }

            // Note: In Rust, (a..a) is the empty iterator.
            node.next
                .iter()
                .skip(node.i + 1)
                .filter_map(|o| o.as_ref())
                .for_each(|n| {
                    self.push(Seed {
                        part,
                        start: 0,
                        node: &**n,
                    });
                });
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L> + std::fmt::Debug> SylowStreamBuilder<S, L, C> {
    /// Returns a new `SylowStreamBuilder`.
    pub fn new() -> SylowStreamBuilder<S, L, C> {
        SylowStreamBuilder {
            mode: flags::NONE,
            tree: FactorNode {
                i: 0,
                rs: [0; L],
                step: 0,
                lim: 0,
                consume: false,
                next: std::array::from_fn(|_| None),
            },
            _phantom: PhantomData,
        }
    }

    /// Adds a flag to the `SylowStreamBuilder`, modifying its yields.
    pub fn add_flag(mut self, mode: u8) -> SylowStreamBuilder<S, L, C> {
        self.mode |= mode;
        self
    }

    /// Adds a target order to this `SylowStreamBuilder`.
    /// The `SylowStream` built from this builder will only yield elements of the orders of
    /// `target`s, or elements of order dividing `target` if `target
    pub fn add_target(mut self, t: [usize; L]) -> SylowStreamBuilder<S, L, C> {
        if t.iter().all(|x| *x == 0) {
            self.mode |= flags::INCLUDE_ONE;
        }

        struct Recur<'a, const L1: usize> {
            f: &'a dyn Fn(&Recur<L1>, &mut FactorNode<L1>, usize, bool),
        }
        let add_target_help = Recur::<L> {
            f: &|recur, node, i, block| {
                if self.mode & flags::LEQ != 0
                    || (t[i] == node.rs[i] && {
                        let mut j = i + 1;
                        loop {
                            if j == L {
                                break true;
                            }
                            if t[j] != 0 {
                                break false;
                            }
                            j += 1;
                        }
                    })
                {
                    node.consume = true;
                }
                let (p, _) = C::FACTORS.prime_powers()[i];
                for j in i..L {
                    let (p2, d2) = C::FACTORS.prime_powers()[j];
                    if t[j] > node.rs[j] {
                        let mut rs = node.rs;
                        rs[j] += 1;
                        let mut lim = C::FACTORS.factor(j);
                        let block = block
                            && ((i == 0 && node.rs[0] == 0)
                                || j == i
                                || (p == 2 && node.rs[node.i] == 1));
                        if block {
                            lim /= 2;
                        }
                        (recur.f)(
                            recur,
                            node.next[j].get_or_insert(Box::new(FactorNode {
                                i: j,
                                rs,
                                step: intpow::<0>(p2, (d2 - 1 - node.rs[j]) as u128),
                                consume: false,
                                lim,
                                next: std::array::from_fn(|_| None),
                            })),
                            j,
                            block,
                        );
                        if self.mode & flags::LEQ == 0 {
                            break;
                        }
                    }
                }
            },
        };
        (add_target_help.f)(
            &add_target_help,
            &mut self.tree,
            0,
            self.mode & flags::NO_UPPER_HALF != 0,
        );

        self
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Iterator for SylowSeqStream<S, L, C> {
    type Item = SylowElem<S, L, C>;

    fn next(&mut self) -> Option<SylowElem<S, L, C>> {
        if let Some(res) = self.buffer.pop() {
            Some(res)
        } else if let Some(top) = self.stack.pop() {
            self.propogate(top, |slf, e| slf.buffer.push(e));
            self.next()
        } else {
            None
        }
    }
}

impl<S, const L: usize, C> SylowParStream<S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L> + Send + Sync,
{
    fn maybe_split(&mut self, stolen: bool) -> Option<Self> {
        if stolen {
            self.splits = rayon::current_num_threads();
        }

        if self.splits == 0 {
            return None;
        }

        let len = self.stack.len();
        if len <= 1 {
            return None;
        }
        let stack = self.stack.split_off(len / 2);
        self.splits /= 2;
        Some(SylowParStream {
            tree: Arc::clone(&self.tree),
            stack,
            splits: self.splits,
            emit_one: false,
        })
    }

    fn work<Con>(&mut self, stolen: bool, consumer: Con) -> Con::Result
    where
        C: Send,
        Con: UnindexedConsumer<SylowElem<S, L, C>>,
    {
        let mut folder = consumer.split_off_left().into_folder();
        if self.emit_one {
            folder = folder.consume(SylowElem::one());
            self.emit_one = false;
        }
        let folder = RefCell::new(Some(folder));

        let mut count = 0;
        loop {
            if let Some(top) = self.stack.pop() {
                self.propogate(top, |_, e| {
                    let mut f = folder.take().unwrap();
                    f = f.consume(e);
                    folder.replace(Some(f));
                });
            } else {
                break;
            }

            if count % 10_000 == 0 {
                let Some(mut split) = self.maybe_split(stolen) else {
                    continue;
                };
                let (r1, r2) = (consumer.to_reducer(), consumer.to_reducer());
                let left_consumer = consumer.split_off_left();

                let (left, right) = rayon::join_context(
                    |ctx| self.work(ctx.migrated(), left_consumer),
                    |ctx| split.work(ctx.migrated(), consumer),
                );
                return r1.reduce(
                    folder.into_inner().unwrap().complete(),
                    r2.reduce(left, right),
                );
            }
            count += 1;
        }
        folder.into_inner().unwrap().complete()
    }
}

impl<S, const L: usize, C> ParallelIterator for SylowParStream<S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L> + Send + Sync,
{
    type Item = SylowElem<S, L, C>;

    fn drive_unindexed<Con>(mut self, consumer: Con) -> Con::Result
    where
        Con: UnindexedConsumer<Self::Item>,
    {
        self.work(false, consumer)
    }
}

impl<S, const L: usize, C> IntoIterator for SylowStreamBuilder<S, L, C>
where
    C: SylowDecomposable<S, L>,
{
    type Item = SylowElem<S, L, C>;
    type IntoIter = SylowSeqStream<S, L, C>;

    fn into_iter(self) -> SylowSeqStream<S, L, C> {
        let mut buffer = Vec::new();
        if (self.mode & flags::INCLUDE_ONE != 0)
            || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0)
        {
            buffer.push(SylowElem::one());
        }
        let mut stream = SylowSeqStream {
            tree: self.tree,
            stack: Vec::new(),
            buffer,
        };
        stream.init_stack(self.mode);
        stream
    }
}

impl<S, const L: usize, C> IntoParallelIterator for SylowStreamBuilder<S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L> + Send + Sync,
{
    type Item = SylowElem<S, L, C>;
    type Iter = SylowParStream<S, L, C>;

    fn into_par_iter(self) -> Self::Iter {
        let mut res = SylowParStream {
            splits: rayon::current_num_threads(),
            stack: Vec::new(),
            emit_one: self.mode & flags::INCLUDE_ONE != 0
                || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0),
            tree: Arc::new(self.tree),
        };
        res.init_stack(self.mode);
        res
    }
}

impl<'a, S, const L: usize, C> SylowStream<'a, S, L, C> for SylowSeqStream<S, L, C>
where
    C: SylowDecomposable<S, L>,
{
    fn push(&mut self, e: Seed<S, L, C>) {
        self.stack.push(e);
    }

    fn stack(&mut self) -> &mut Vec<Seed<S, L, C>> {
        &mut self.stack
    }

    fn tree(&self) -> &FactorNode<L> {
        &self.tree
    }
}

impl<'a, S, const L: usize, C> SylowStream<'a, S, L, C> for SylowParStream<S, L, C>
where
    S: Send + Sync,
    C: SylowDecomposable<S, L>,
{
    fn push(&mut self, e: Seed<S, L, C>) {
        self.stack.push(e);
    }

    fn stack(&mut self) -> &mut Vec<Seed<S, L, C>> {
        &mut self.stack
    }

    fn tree(&self) -> &FactorNode<L> {
        &*self.tree
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Clone for Seed<S, L, C> {
    fn clone(&self) -> Seed<S, L, C> {
        Seed {
            part: self.part,
            start: self.start,
            node: self.node,
        }
    }
}
impl<S, const L: usize, C: SylowDecomposable<S, L>> Copy for Seed<S, L, C> {}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Clone for SylowStreamBuilder<S, L, C> {
    fn clone(&self) -> Self {
        SylowStreamBuilder {
            mode: self.mode,
            tree: self.tree.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Clone for SylowSeqStream<S, L, C> {
    fn clone(&self) -> SylowSeqStream<S, L, C> {
        SylowSeqStream {
            stack: self.stack.clone(),
            buffer: self.buffer.clone(),
            tree: self.tree.clone(),
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S, L>> Clone for SylowParStream<S, L, C>
where
    S: Send + Sync,
{
    fn clone(&self) -> SylowParStream<S, L, C> {
        SylowParStream {
            stack: self.stack.clone(),
            splits: self.splits,
            emit_one: self.emit_one,
            tree: Arc::clone(&self.tree),
        }
    }
}

unsafe impl<S, const L: usize, C: SylowDecomposable<S, L>> Send for Seed<S, L, C> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::fp::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    impl Factor<Phantom, 2> for FpNum<7> {
        const FACTORS: Factorization<2> = Factorization::new(&[(2, 1), (3, 1)]);
    }

    impl Factor<Phantom, 3> for FpNum<61> {
        const FACTORS: Factorization<3> = Factorization::new(&[(2, 2), (3, 1), (5, 1)]);
    }

    impl Factor<Phantom, 7> for FpNum<BIG_P> {
        const FACTORS: Factorization<7> = Factorization::new(&[
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
    }

    impl Factor<Phantom, 3> for FpNum<271> {
        const FACTORS: Factorization<3> = Factorization::new(&[(2, 1), (3, 3), (5, 1)]);
    }

    impl Factor<Phantom, 4> for FpNum<13928643> {
        const FACTORS: Factorization<4> = Factorization::new(&[(2, 1), (7, 2), (13, 2), (29, 2)]);
    }

    impl Factor<Phantom, 2> for FpNum<41> {
        const FACTORS: Factorization<2> = Factorization::new(&[(2, 3), (5, 1)]);
    }

    #[test]
    pub fn test_make_stream_seq() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new()
            .add_target([1, 0])
            .into_iter()
            .filter_map(|s| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_seq() {
        let stream = SylowStreamBuilder::new().add_target([1, 0, 0]).into_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<61>>> = stream.collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(!x.is_one());
        x = x.pow(2);
        assert!(x.is_one());

        let mut count = 0;
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
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

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target([0, 1, 0])
            .into_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_leq_seq() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_target([2, 1, 0])
            .into_iter()
            .count();
        assert_eq!(count, 12);
    }

    #[test]
    pub fn test_generates_big_seq() {
        let stream = SylowStreamBuilder::new()
            .add_target([0, 0, 0, 2, 0, 0, 0])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
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
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new().add_target([0, 2, 1]);
        let stream_all = builder.into_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([0, 2, 1])
            .into_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_target([1, 0, 0])
            .add_target([0, 1, 0])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::new()
            .add_flag(flags::LEQ)
            .add_target([1, 1, 0])
            .add_target([0, 2, 0])
            .add_target([0, 0, 1])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_seq() {
        let coords: Vec<SylowElem<Phantom, 4, FpNum<13928643>>> = SylowStreamBuilder::new()
            .add_flag(flags::LEQ)
            .add_target([0, 1, 1, 0])
            .into_iter()
            .collect();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_seq() {
        let mut count = 0;
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_target([2, 0, 1])
            .into_iter()
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.multiply(&x);
                assert!(!x.is_one());
                count += 1;
            });
        assert_eq!(count, 18);
    }

    #[test]
    pub fn test_make_stream_par() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new()
            .add_target([1, 0])
            .into_par_iter()
            .filter_map(|s| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_par() {
        let stream = SylowStreamBuilder::new()
            .add_target([1, 0, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<61>>> = stream.collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(!x.is_one());
        x = x.pow(2);
        assert!(x.is_one());

        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
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

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target([0, 1, 0])
            .into_par_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_generates_big_par() {
        let stream = SylowStreamBuilder::new()
            .add_target([0, 0, 0, 2, 0, 0, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
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
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new().add_target([0, 2, 1]);
        let stream_all = builder.into_par_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([0, 2, 1])
            .into_par_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_par() {
        let stream = SylowStreamBuilder::new()
            .add_target([1, 0, 0])
            .add_target([0, 1, 0])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::new()
            .add_flag(flags::LEQ)
            .add_target([1, 1, 0])
            .add_target([0, 2, 0])
            .add_target([0, 0, 1])
            .into_par_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<271>>> = stream.collect();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_par() {
        let coords: Vec<SylowElem<Phantom, 4, FpNum<13928643>>> = SylowStreamBuilder::new()
            .add_flag(flags::LEQ)
            .add_target([0, 1, 1, 0])
            .into_par_iter()
            .collect();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_par() {
        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_target([2, 0, 1])
            .into_par_iter()
            .for_each(|mut x| {
                assert!(!x.is_one());
                x = x.multiply(&x);
                assert!(!x.is_one());
                count.fetch_add(1, Ordering::Relaxed);
            });
        assert_eq!(count.into_inner(), 18);
    }

    #[test]
    pub fn test_no_parabolic_no_upper_half_seq() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([2, 0, 1])
            .into_iter()
            .count();
        assert_eq!(count, 9);
        // assert!(false);
    }

    #[test]
    pub fn test_no_parabolic_no_upper_half_par() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([2, 0, 1])
            .into_par_iter()
            .count();
        assert_eq!(count, 9);
    }

    #[test]
    pub fn test_subordinate_target() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target([0, 1, 0])
            .add_target([0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 10);
    }

    #[test]
    pub fn test_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 4);
    }

    #[test]
    pub fn test_propogates_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([3, 1])
            .into_iter()
            .count();
        assert_eq!(count, 8);

        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target([1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 2);
    }
}
