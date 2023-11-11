use rayon::iter::plumbing::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use std::cell::RefCell;
use std::marker::{PhantomData};
use std::sync::Arc;

use crate::numbers::*;
use crate::streams::FactorStream;
use libbgs_util::*;

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
pub struct SylowStreamBuilder<S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug> {
    mode: u8,
    tree: Box<FactorNode<L>>,
    _phantom: PhantomData<(S, C)>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements in parallel on multiple threads.
pub struct SylowParStream<
    S: Send + Sync,
    const L: usize,
    C: SylowDecomposable<S> + std::fmt::Debug,
> {
    stack: Vec<Seed<S, L, C>>,
    splits: usize,
    buffer: Vec<Output<S, L, C>>,
    tree: Arc<FactorNode<L>>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements sequentially on a single thread.
pub struct SylowSeqStream<S, const L: usize, C: SylowDecomposable<S>> {
    stack: Vec<Seed<S, L, C>>,
    buffer: Vec<Output<S, L, C>>,
    tree: Arc<FactorNode<L>>,
}

#[derive(Debug)]
struct Seed<S, const L: usize, C: SylowDecomposable<S>> {
    part: SylowElem<S, L, C>,
    start: u128,
    node: *const FactorNode<L>,
}

#[derive(Clone, Debug)]
struct GenData {
    consume: bool,
    step: u128,
    lim: u128,
}

type FactorNode<const L: usize> = FactorTrie<L, GenData>;
type Output<S, const L: usize, C> = (SylowElem<S, L, C>, [usize; L]);

trait SylowStream<'a, S, const L: usize, C>
where
    Self: Sized,
    C: SylowDecomposable<S>,
{
    fn push(&mut self, e: Seed<S, L, C>);
    fn stack(&mut self) -> &mut Vec<Seed<S, L, C>>;
    fn tree(&self) -> &FactorNode<L>;

    fn init_stack(&mut self, mode: u8) {
        for i in 0..L {
            let Some(n) = self.tree().child(i) else {
                continue;
            };

            let seed = Seed {
                part: SylowElem::ONE,
                start: 0,
                node: &*n,
            };

            let (p, _) = C::FACTORS[i];
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
        F: FnMut(&mut Self, Output<S, L, C>),
    {
        let node = unsafe { &*seed.node };
        let (p, _) = C::FACTORS[node.index()];

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
            let tmp = seed.part.coords[node.index()] + j * node.data.step;
            if tmp > node.data.lim {
                break;
            }
            let mut part = seed.part;
            part.coords[node.index()] = tmp;

            if let Some(n) = node.child(node.index()) {
                self.push(Seed {
                    part,
                    start: 0,
                    node: &*n,
                });
            }

            // Next, create new seeds by moving to the next prime power,
            // but only if we are *done* with this prime power.
            if j == 0 {
                continue;
            }
            if node.data.consume {
                consume(self, (part, *node.ds()));
            }

            node.children()
                .iter()
                .skip(node.index() + 1)
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

impl<S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug> SylowStreamBuilder<S, L, C> {
    /// Returns a new `SylowStreamBuilder`.
    pub fn new() -> SylowStreamBuilder<S, L, C> {
        SylowStreamBuilder {
            mode: flags::NONE,
            tree: Box::new(FactorNode::new(GenData {
                step: 0,
                lim: 0,
                consume: false,
            })),
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
    pub fn add_target(mut self, t: &[usize; L]) -> SylowStreamBuilder<S, L, C> {
        if t.iter().all(|x| *x == 0) {
            self.mode |= flags::INCLUDE_ONE;
        }

        struct Adder<'a, S1, const L1: usize, C1> {
            block: bool,
            mode: u8,
            t: &'a [usize; L1],
            _phantom: PhantomData<(S1, C1)>,
        }

        impl<'a, S1, const L1: usize, C1> MutFactorVisitor<L1, GenData> for Adder<'a, S1, L1, C1> 
        where
            C1: Factor<S1>,
        {
            fn visit_mut(&mut self, node: &mut FactorTrie<L1, GenData>) {
                node.data.consume |= self.mode & flags::LEQ != 0
                    || (self.t[node.index()] == node.ds()[node.index()] && {
                        let mut j = node.index() + 1;
                        loop {
                            if j == L1 {
                                break true;
                            }
                            if self.t[j] != 0 {
                                break false;
                            }
                            j += 1;
                        }
                    });
                let (p, _) = C1::FACTORS.prime_powers()[node.index()];
                for j in node.index()..L1 {
                    let (p2, d2) = C1::FACTORS.prime_powers()[j];
                    if self.t[j] > node.ds()[j] {
                        let mut rs = node.ds().clone();
                        rs[j] += 1;
                        let mut lim = C1::FACTORS.factor(j);
                        let block = self.block
                            && ((node.index() == 0 && node.ds()[0] == 0)
                                || j == node.index()
                                || (p == 2 && node.ds()[node.index()] == 1));
                        if block {
                            lim /= 2;
                        }
                        let child = node.get_or_new_child(j, GenData {
                            step: intpow::<0>(p2, (d2 - 1 - node.ds()[j]) as u128),
                            consume: false,
                            lim,
                        });
                        Adder {
                            block,
                            ..*self
                        }.visit_mut(child);
                        if self.mode & flags::LEQ == 0 {
                            break;
                        }
                    }
                }
            }
        }

        Adder {
            block: self.mode & flags::NO_UPPER_HALF != 0,
            mode: self.mode,
            t,
            _phantom: PhantomData::<(S, C)>,
        }.visit_mut(&mut self.tree);

        self
    }

    pub fn add_quotient(mut self, q: [usize; L]) -> Self {
        struct Quotienter<const L1: usize> {
            lims: [u128; L1],
        }

        impl<const L1: usize> Quotienter<L1> {
            fn visit_mut(&self, node: &mut FactorTrie<L1, GenData>) {
                if node.data.step > self.lims[node.index()] {
                    node.data.lim = self.lims[node.index()];
                }
                for j in node.index()..L1 {
                    if let Some(ref mut child) = node.child_mut(j) {
                        self.visit_mut(child);
                    }
                }
            }
        }

        Quotienter {
            lims: std::array::from_fn(|i| {
                let (p, d) = C::FACTORS.prime_powers()[i];
                if q[i] <= d {
                    intpow::<0>(p, (d - q[i]) as u128) - 1
                } else { 0 }
            })
        }.visit_mut(&mut self.tree);

        self
    }

    /// Add all the targets yielded by this `FactorStream`.
    pub fn add_targets_from_factors(self, stream: FactorStream) -> Self {
        stream.map(|v| v.try_into().unwrap())
            .fold(self, |b, x| b.add_target(&x))
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> SylowSeqStream<S, L, C> {
    /// Converts a sequential Sylow stream into a parallel one.
    pub fn parallelize(self) -> SylowParStream<S, L, C>
    where
        S: Send + Sync,
    {
        SylowParStream {
            stack: self.stack,
            buffer: self.buffer,
            splits: rayon::current_num_threads(),
            tree: Arc::from(self.tree),
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> Iterator for SylowSeqStream<S, L, C> {
    type Item = (SylowElem<S, L, C>, [usize; L]);

    fn next(&mut self) -> Option<(SylowElem<S, L, C>, [usize; L])> {
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
    C: SylowDecomposable<S> + Send + Sync,
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
            buffer: Vec::new(),
        })
    }

    fn work<Con>(&mut self, stolen: bool, consumer: Con) -> Con::Result
    where
        C: Send,
        Con: UnindexedConsumer<Output<S, L, C>>,
    {
        let mut folder = consumer.split_off_left().into_folder();
        while let Some(buf) = self.buffer.pop() {
            folder = folder.consume(buf);
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
    C: SylowDecomposable<S> + Send + Sync,
{
    type Item = Output<S, L, C>;

    fn drive_unindexed<Con>(mut self, consumer: Con) -> Con::Result
    where
        Con: UnindexedConsumer<Self::Item>,
    {
        self.work(false, consumer)
    }
}

impl<S, const L: usize, C> IntoIterator for SylowStreamBuilder<S, L, C>
where
    C: SylowDecomposable<S>,
{
    type Item = Output<S, L, C>;
    type IntoIter = SylowSeqStream<S, L, C>;

    fn into_iter(self) -> SylowSeqStream<S, L, C> {
        let mut buffer = Vec::new();
        if (self.mode & flags::INCLUDE_ONE != 0)
            || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0)
        {
            buffer.push((SylowElem::ONE, [0; L]));
        }
        let mut stream = SylowSeqStream {
            tree: Arc::from(self.tree),
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
    C: SylowDecomposable<S> + Send + Sync,
{
    type Item = Output<S, L, C>;
    type Iter = SylowParStream<S, L, C>;

    fn into_par_iter(self) -> Self::Iter {
        let mut res = SylowParStream {
            splits: rayon::current_num_threads(),
            stack: Vec::new(),
            buffer: if self.mode & flags::INCLUDE_ONE != 0
                || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0) {
                    vec![(SylowElem::ONE, [0; L])]
                } else {
                    Vec::new()
                },
            tree: Arc::from(self.tree),
        };
        res.init_stack(self.mode);
        res
    }
}

impl<'a, S, const L: usize, C> SylowStream<'a, S, L, C> for SylowSeqStream<S, L, C>
where
    C: SylowDecomposable<S>,
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
    C: SylowDecomposable<S>,
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

impl<S, const L: usize, C: SylowDecomposable<S>> Clone for Seed<S, L, C> {
    fn clone(&self) -> Seed<S, L, C> {
        Seed {
            part: self.part,
            start: self.start,
            node: self.node,
        }
    }
}
impl<S, const L: usize, C: SylowDecomposable<S>> Copy for Seed<S, L, C> {}

impl<S, const L: usize, C: SylowDecomposable<S>> Clone for SylowStreamBuilder<S, L, C> {
    fn clone(&self) -> Self {
        SylowStreamBuilder {
            mode: self.mode,
            tree: self.tree.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> Clone for SylowSeqStream<S, L, C> {
    fn clone(&self) -> SylowSeqStream<S, L, C> {
        SylowSeqStream {
            stack: self.stack.clone(),
            buffer: self.buffer.clone(),
            tree: self.tree.clone(),
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> Clone for SylowParStream<S, L, C>
where
    S: Send + Sync,
{
    fn clone(&self) -> SylowParStream<S, L, C> {
        SylowParStream {
            stack: self.stack.clone(),
            splits: self.splits,
            buffer: self.buffer.clone(),
            tree: Arc::clone(&self.tree),
        }
    }
}

unsafe impl<S, const L: usize, C: SylowDecomposable<S>> Send for Seed<S, L, C> {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq)]
    struct Phantom {}

    impl Factor<Phantom> for FpNum<7> {
        const FACTORS: Factorization = Factorization::new(&[(2, 1), (3, 1)]);
    }

    impl Factor<Phantom> for FpNum<61> {
        const FACTORS: Factorization = Factorization::new(&[(2, 2), (3, 1), (5, 1)]);
    }

    impl Factor<Phantom> for FpNum<BIG_P> {
        const FACTORS: Factorization = Factorization::new(&[
            (2, 1),
            (7, 1),
            (13, 1),
            (29, 2),
            (43, 1),
            (705737, 1),
            (215288719, 1),
        ]);
    }

    impl Factor<Phantom> for FpNum<271> {
        const FACTORS: Factorization = Factorization::new(&[(2, 1), (3, 3), (5, 1)]);
    }

    impl Factor<Phantom> for FpNum<13928643> {
        const FACTORS: Factorization = Factorization::new(&[(2, 1), (7, 2), (13, 2), (29, 2)]);
    }

    impl Factor<Phantom> for FpNum<41> {
        const FACTORS: Factorization = Factorization::new(&[(2, 3), (5, 1)]);
    }

    #[test]
    pub fn test_make_stream_seq() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new()
            .add_target(&[1, 0])
            .into_iter()
            .filter_map(|(s, _)| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_seq() {
        let stream = SylowStreamBuilder::new().add_target(&[1, 0, 0]).into_iter();
        let coords: Vec<SylowElem<Phantom, 3, FpNum<61>>> = stream.map(|(a, _)| a).collect();
        assert_eq!(coords.len(), 1);
        let mut x = coords[0].clone();
        assert!(x != SylowElem::ONE);
        x = x.pow(2);
        assert!(x == SylowElem::ONE);

        let mut count = 0;
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[2, 0, 0])
            .into_iter()
            .for_each(|(mut x, _)| {
                count += 1;
                for _ in 1..3 {
                    let y = x.clone();
                    assert!(x != SylowElem::ONE);
                    x = x.multiply(&y);
                }
                assert!(x == SylowElem::ONE);
            });
        assert_eq!(count, 2);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[0, 1, 0])
            .into_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_leq_seq() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_target(&[2, 1, 0])
            .into_iter()
            .count();
        assert_eq!(count, 12);
    }

    #[test]
    pub fn test_generates_big_seq() {
        let stream = SylowStreamBuilder::new()
            .add_target(&[0, 0, 0, 2, 0, 0, 0])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.map(|(a, _)| a).collect();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
            .add_target(&[0, 0, 0, 0, 0, 1, 0])
            .into_iter()
            .take(2)
            .for_each(|(mut x, _)| {
                assert!(x != SylowElem::ONE);
                x = x.pow(705737);
                assert!(x == SylowElem::ONE);
            });
    }

    #[test]
    pub fn test_generates_medium_seq() {
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new().add_target(&[0, 2, 1]);
        let stream_all = builder.into_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 2, 1])
            .into_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_target(&[1, 0, 0])
            .add_target(&[0, 1, 0])
            .into_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::LEQ)
            .add_target(&[1, 1, 0])
            .add_target(&[0, 2, 0])
            .add_target(&[0, 0, 1])
            .into_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_seq() {
        let coords = SylowStreamBuilder::<Phantom, 4, FpNum<13928643>>::new()
            .add_flag(flags::LEQ)
            .add_target(&[0, 1, 1, 0])
            .into_iter()
            .collect::<Vec<_>>();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_seq() {
        let mut count = 0;
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_target(&[2, 0, 1])
            .into_iter()
            .for_each(|(mut x, _)| {
                assert!(x != SylowElem::ONE);
                x = x.multiply(&x);
                assert!(x != SylowElem::ONE);
                count += 1;
            });
        assert_eq!(count, 18);
    }

    #[test]
    pub fn test_make_stream_par() {
        let g = SylowDecomp::<Phantom, 2, FpNum<7>>::new();
        let res: Vec<FpNum<7>> = SylowStreamBuilder::new()
            .add_target(&[1, 0])
            .into_par_iter()
            .filter_map(|(s, _)| s.to_product(&g).into())
            .collect();
        assert_eq!(res, vec![6]);
    }

    #[test]
    pub fn test_generates_small_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[1, 0, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 1);
        let (mut x, _) = coords[0].clone();
        assert!(x != SylowElem::ONE);
        x = x.pow(2);
        assert!(x == SylowElem::ONE);

        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[2, 0, 0])
            .into_par_iter()
            .for_each(|(mut x, _)| {
                count.fetch_add(1, Ordering::Relaxed);
                for _ in 1..3 {
                    let y = x.clone();
                    assert!(x != SylowElem::ONE);
                    x = x.multiply(&y);
                }
                assert!(x == SylowElem::ONE);
            });
        assert_eq!(count.into_inner(), 2);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[0, 1, 0])
            .into_par_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_generates_big_par() {
        let stream = SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
            .add_target(&[0, 0, 0, 2, 0, 0, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
            .add_target(&[0, 0, 0, 0, 0, 1, 0])
            .into_par_iter()
            .take_any(2)
            .for_each(|(mut x, _)| {
                assert!(x != SylowElem::ONE);
                x = x.pow(705737);
                assert!(x == SylowElem::ONE);
            });
    }

    #[test]
    pub fn test_generates_medium_par() {
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new().add_target(&[0, 2, 1]);
        let stream_all = builder.into_par_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 2, 1])
            .into_par_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_target(&[1, 0, 0])
            .add_target(&[0, 1, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>>::new()
            .add_flag(flags::LEQ)
            .add_target(&[1, 1, 0])
            .add_target(&[0, 2, 0])
            .add_target(&[0, 0, 1])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 16);
    }

    #[test]
    pub fn test_multiple_targets_2_par() {
        let coords = SylowStreamBuilder::<Phantom, 4, FpNum<13928643>>::new()
            .add_flag(flags::LEQ)
            .add_target(&[0, 1, 1, 0])
            .into_par_iter()
            .collect::<Vec<_>>();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_par() {
        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_target(&[2, 0, 1])
            .into_par_iter()
            .for_each(|(mut x, _)| {
                assert!(x != SylowElem::ONE);
                x = x.multiply(&x);
                assert!(x != SylowElem::ONE);
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
            .add_target(&[2, 0, 1])
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
            .add_target(&[2, 0, 1])
            .into_par_iter()
            .count();
        assert_eq!(count, 9);
    }

    #[test]
    pub fn test_subordinate_target() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_target(&[0, 1, 0])
            .add_target(&[0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 10);
    }

    #[test]
    pub fn test_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 4);
    }

    #[test]
    pub fn test_propogates_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[3, 1])
            .into_iter()
            .count();
        assert_eq!(count, 8);

        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 2);
    }
}
