use rayon::iter::plumbing::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use std::cell::RefCell;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::numbers::*;
use crate::streams::DivisorStream;
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
pub struct SylowStreamBuilder<S, const L: usize, C: SylowDecomposable<S>, T> {
    mode: u8,
    tree: Box<FactorTrie<S, L, C, (Consume, T)>>,
    quotient: Option<[usize; L]>,
    _phantom: PhantomData<(S, C)>,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements in parallel on multiple threads.
pub struct SylowParStream<S: Send + Sync, const L: usize, C: SylowDecomposable<S>, T>
{
    stream: SylowStream<S, L, C, T>,
    splits: usize,
}

/// A stream yielding elements of particular orders, as their Sylow decompositions.
/// Generates the elements sequentially on a single thread.
pub struct SylowStream<S, const L: usize, C: SylowDecomposable<S>, T> {
    stack: Vec<Seed<S, L, C, T>>,
    buffer: Vec<Output<S, L, C, T>>,
    tree: Arc<FactorTrie<S, L, C, (GenData, T)>>,
}

#[derive(Debug)]
struct Seed<S, const L: usize, C: SylowDecomposable<S>, T> {
    part: SylowElem<S, L, C>,
    start: u128,
    node: *const FactorTrie<S, L, C, (GenData, T)>,
}

#[derive(Clone, Debug)]
struct GenData {
    consume: Consume,
    step: u128,
    lim: u128,
}

type Output<S, const L: usize, C, T> = (SylowElem<S, L, C>, T);
#[derive(Clone, Debug)]
struct Consume {
    this: bool,
    descendants: usize,
}
impl Default for Consume {
    fn default() -> Self {
        Consume { this: false, descendants: 0 }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug> SylowStreamBuilder<S, L, C, ()> {
    /// Returns a new `SylowStreamBuilder`.
    pub fn new() -> SylowStreamBuilder<S, L, C, ()> {
        SylowStreamBuilder {
            mode: flags::NONE,
            tree: Box::new(FactorTrie::new().map(&|_: (), _, _| (Consume::default(), ()))),
            quotient: None,
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug> SylowStreamBuilder<S, L, C, [u128; L]> {
    /// Returns a new `SylowStreamBuilder`, which will return both elements and their orders.
    pub fn new_with_orders() -> SylowStreamBuilder<S, L, C, [usize; L]> {
        SylowStreamBuilder {
            mode: flags::NONE,
            tree: Box::new(FactorTrie::<S, L, C, ()>::new().map(&|_, ds, _| (Consume::default(), *ds))),
            quotient: None,
            _phantom: PhantomData,
        }
    }
}

impl<'a, S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug, T> SylowStreamBuilder<S, L, C, &'a T> {
    /// Creates a new `SylowStreamBuilder` with a "parallel" trie to that given here.
    pub fn new_with_trie(trie: &'a FactorTrie<S, L, C, T>) -> SylowStreamBuilder<S, L, C, &'a T> {
        SylowStreamBuilder {
            mode: flags::NONE,
            tree: Box::new(trie.as_ref().map(&|t, _, _| (Consume::default(), t))),
            quotient: None,
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S> + std::fmt::Debug, T> SylowStreamBuilder<S, L, C, T> {
    /// Adds a flag to the `SylowStreamBuilder`, modifying its yields.
    pub fn add_flag(mut self, mode: u8) -> SylowStreamBuilder<S, L, C, T> {
        self.mode |= mode;
        self
    }

    /// Adds a target order to this `SylowStreamBuilder`.
    /// The `SylowStream` built from this builder will only yield elements of the orders of
    /// `target`s, or elements of order dividing `target` if `target
    pub fn add_target(mut self, t: &[usize; L]) -> SylowStreamBuilder<S, L, C, T> {
        if t.iter().all(|x| *x == 0) {
            self.mode |= flags::INCLUDE_ONE;
        }

        struct Adder<'a, const L1: usize> {
            mode: u8,
            t: &'a [usize; L1],
        }

        impl<'a, const L1: usize> Adder<'a, L1> {
            fn visit_mut<S1, C1, T1>(&mut self, node: &mut FactorTrie<S1, L1, C1, (Consume, T1)>) -> usize {
                node.data.0.this |= self.mode & flags::LEQ != 0
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
                for j in node.index()..L1 {
                    if self.t[j] > node.ds()[j] {
                        let Some(child) = node.child_mut(j) else { 
                            panic!("Tried to add a target which does not exist in this trie!"); 
                        };
                        node.data.0.descendants = self.visit_mut(child);
                        if self.mode & flags::LEQ == 0 {
                            break;
                        }
                    }
                }
                node.data.0.descendants + if node.data.0.this { 1 } else { 0 }
            }
        }

        Adder { mode: self.mode, t }.visit_mut(&mut self.tree);

        self
    }

    /// Adds the maximal divisors beneath `limit` to the `FactorTrie` and sets the `LEQ` flag.
    pub fn add_targets_leq(self, limit: u128) -> Self {
        DivisorStream::new(C::FACTORS.factors(), limit, true)
            .into_iter()
            .map(|v| v.try_into().unwrap())
            .fold(self, |b, x| b.add_target(&x))
            .add_flag(flags::LEQ)
    }

    /// Remove the target, so elements of that order will not be generated.
    pub fn remove_target(mut self, t: &[usize; L]) -> Self {
        if t.iter().all(|x| *x == 0) {
            self.mode |= flags::INCLUDE_ONE;
        }

        struct Remover<'a, const L1: usize>(&'a [usize; L1]);
        impl<'a, const L1: usize> Remover<'a, L1> {
            fn visit_mut<S1, C1, T1>(&mut self, node: &mut FactorTrie<S1, L1, C1, (Consume, T1)>) -> bool {
                for j in node.index()..L1 {
                    if self.0[j] > node.ds()[j] {
                        let Some(child) = node.child_mut(j) else { 
                            panic!("Could not find child while removing target."); 
                        };
                        if self.visit_mut(child) {
                            node.data.0.descendants -= 1;
                            return true;
                        }
                        return false;
                    }
                }
                node.data.0.this = false;
                true
            }
        }
        Remover(t).visit_mut(&mut self.tree);
        self
    }

    /// Guarantees that this stream will only ever yield one representative of the cosets of the
    /// quotient.
    pub fn set_quotient(mut self, q: Option<[usize; L]>) -> Self {
        self.quotient = q;
        self
    }

    /// Add all the targets yielded by this `DivisorStream`.
    pub fn add_targets_from_factors(self, stream: DivisorStream) -> Self {
        stream
            .map(|v| v.try_into().unwrap())
            .fold(self, |b, x| b.add_target(&x))
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>, T> SylowStream<S, L, C, T> {
    /// Converts a sequential Sylow stream into a parallel one.
    pub fn parallelize(self) -> SylowParStream<S, L, C, T>
    where
        S: Send + Sync,
    {
        SylowParStream {
            stream: self,
            splits: rayon::current_num_threads(),
        }
    }

    fn propagate<F>(&mut self, seed: Seed<S, L, C, T>, mut consume: F)
    where
        Self: Sized,
        T: Clone,
        F: FnMut(&mut Self, Output<S, L, C, T>),
    {
        let node = unsafe { &*seed.node };
        let (p, _) = C::FACTORS[node.index()];

        // First, create new seeds by incrementing
        // the current power.
        let mut stop = p;

        if stop - seed.start > STACK_ADDITION_LIMIT as u128 {
            self.stack.push(Seed {
                start: seed.start + STACK_ADDITION_LIMIT as u128,
                ..seed
            });
            stop = seed.start + STACK_ADDITION_LIMIT as u128;
        }

        for j in seed.start..stop {
            let tmp = seed.part.coords[node.index()] + j * node.data.0.step;
            if tmp > node.data.0.lim {
                break;
            }
            let mut part = seed.part;
            part.coords[node.index()] = tmp;

            if let Some(n) = node.child(node.index()) {
                self.stack.push(Seed {
                    part,
                    start: 0,
                    node: n,
                });
            }

            // Next, create new seeds by moving to the next prime power,
            // but only if we are *done* with this prime power.
            if j == 0 {
                continue;
            }
            if node.data.0.consume.this {
                consume(self, (part, node.data.1.clone()));
            }

            for i in (node.index() + 1)..L {
                let Some(n) = node.child(i) else { continue; };
                if n.data.0.consume.this || n.data.0.consume.descendants >= 1 {
                    self.stack.push(Seed {
                        part,
                        start: 0,
                        node: &*n,
                    });
                }
            }
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>, T> Iterator for SylowStream<S, L, C, T> 
where
    T: Clone,
{
    type Item = (SylowElem<S, L, C>, T);

    fn next(&mut self) -> Option<(SylowElem<S, L, C>, T)> {
        if let Some(res) = self.buffer.pop() {
            Some(res)
        } else if let Some(top) = self.stack.pop() {
            self.propagate(top, |slf, e| slf.buffer.push(e));
            self.next()
        } else {
            None
        }
    }
}

impl<S, const L: usize, C, T> SylowParStream<S, L, C, T>
where
    S: Send + Sync,
    C: SylowDecomposable<S> + Send + Sync,
    T: Clone + Send + Sync,
{
    fn maybe_split(&mut self, stolen: bool) -> Option<Self> {
        if stolen {
            self.splits = rayon::current_num_threads();
        }

        if self.splits == 0 {
            return None;
        }

        let len = self.stream.stack.len();
        if len <= 1 {
            return None;
        }
        let stack = self.stream.stack.split_off(len / 2);
        self.splits /= 2;
        Some(SylowParStream {
            stream: SylowStream {
                tree: Arc::clone(&self.stream.tree),
                stack,
                buffer: Vec::new(),
            },
            splits: self.splits,
        })
    }

    fn work<Con>(&mut self, stolen: bool, consumer: Con) -> Con::Result
    where
        Con: UnindexedConsumer<Output<S, L, C, T>>,
    {
        let mut folder = consumer.split_off_left().into_folder();
        while let Some(buf) = self.stream.buffer.pop() {
            folder = folder.consume(buf);
        }
        let folder = RefCell::new(Some(folder));

        let mut count = 0;
        while let Some(top) = self.stream.next() {
            let mut f = folder.take().unwrap();
            f = f.consume(top);
            folder.replace(Some(f));

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

impl<S, const L: usize, C, T> ParallelIterator for SylowParStream<S, L, C, T>
where
    S: Send + Sync,
    C: SylowDecomposable<S> + Send + Sync,
    T: Clone + Send + Sync,
{
    type Item = Output<S, L, C, T>;

    fn drive_unindexed<Con>(mut self, consumer: Con) -> Con::Result
    where
        Con: UnindexedConsumer<Self::Item>,
    {
        self.work(false, consumer)
    }
}

impl<S, const L: usize, C, T> IntoIterator for SylowStreamBuilder<S, L, C, T>
where
    C: SylowDecomposable<S>,
    T: Clone,
{
    type Item = Output<S, L, C, T>;
    type IntoIter = SylowStream<S, L, C, T>;

    fn into_iter(self) -> SylowStream<S, L, C, T> {
        let mut tree = self.tree.map(&|consume, ds: &[usize; L], i| {
            let (p, d) = C::FACTORS[i];
            (GenData {
                consume: consume.0,
                step: intpow::<0>(p, (d - ds[i]) as u128),
                lim: 0,
            }, consume.1)
        });
        struct Limiter<const L1: usize> {
            block: bool,
            lims: [u128; L1],
        }
        impl<const L1: usize> Limiter<L1> {
            fn visit_mut<S1, C1, T1>(&mut self, node: &mut FactorTrie<S1, L1, C1, (GenData, T1)>) 
            where
                C1: Factor<S1>,
            {
                let (p, _) = C1::FACTORS[node.index()];
                node.data.0.lim = self.lims[node.index()];
                if self.block {
                    node.data.0.lim /= 2;
                }
                for j in node.index()..L1 {
                    Limiter {
                        block: ((p == 2 && node.ds()[0] <= 1) || j == node.index()) && self.block,
                        ..*self
                    }
                    .visit_mut({
                        let Some(child) = node.child_mut(j) else {
                            continue;
                        };
                        child
                    });
                }
            }
        }
        let q = self.quotient.unwrap_or([0; L]);
        Limiter {
            block: self.mode & flags::NO_UPPER_HALF != 0,
            lims: std::array::from_fn(|i| {
                let (p, d) = C::FACTORS[i];
                if q[i] <= d {
                    intpow::<0>(p, (d - q[i]) as u128)
                } else {
                    0
                }
            }),
        }
        .visit_mut(&mut tree);
        let mut stream = SylowStream {
            stack: Vec::new(),
            buffer: if (self.mode & flags::INCLUDE_ONE != 0)
                || (self.mode & flags::LEQ != 0 && self.mode & flags::NO_PARABOLIC == 0)
            {
                vec![(SylowElem::ONE, tree.data.1.clone())]
            } else {
                Vec::<(SylowElem<_, L, _>, _)>::new()
            },
            tree: Arc::from(tree),
        };

        for i in 0..L {
            let Some(n) = stream.tree.child(i) else {
                continue;
            };

            if !n.data.0.consume.this && n.data.0.consume.descendants == 0 {
                continue;
            }

            let seed = Seed {
                part: SylowElem::ONE,
                start: 0,
                node: n,
            };

            let (p, _) = C::FACTORS[i];
            if self.mode & flags::NO_PARABOLIC != 0 && p == 2 {
                stream.propagate(seed, |_, _| {});
            } else {
                stream.stack.push(seed);
            }
        }
        stream
    }
}

impl<S, const L: usize, C, T> IntoParallelIterator for SylowStreamBuilder<S, L, C, T>
where
    S: Send + Sync,
    C: SylowDecomposable<S> + Send + Sync,
    T: Clone + Send + Sync,
{
    type Item = Output<S, L, C, T>;
    type Iter = SylowParStream<S, L, C, T>;

    fn into_par_iter(self) -> Self::Iter {
        SylowParStream {
            stream: self.into_iter(),
            splits: rayon::current_num_threads(),
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>, T: Clone> Clone for Seed<S, L, C, T> {
    fn clone(&self) -> Seed<S, L, C, T> {
        Seed { ..*self }
    }
}
impl<S, const L: usize, C: SylowDecomposable<S>, T: Copy> Copy for Seed<S, L, C, T> {}

impl<S, const L: usize, C: SylowDecomposable<S>, T: Clone> Clone for SylowStreamBuilder<S, L, C, T> {
    fn clone(&self) -> Self {
        SylowStreamBuilder {
            tree: self.tree.clone(),
            ..*self
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>, T: Clone> Clone for SylowStream<S, L, C, T> {
    fn clone(&self) -> SylowStream<S, L, C, T> {
        SylowStream {
            stack: self.stack.clone(),
            buffer: self.buffer.clone(),
            tree: self.tree.clone(),
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>, T: Clone> Clone for SylowParStream<S, L, C, T>
where
    S: Send + Sync,
{
    fn clone(&self) -> SylowParStream<S, L, C, T> {
        SylowParStream {
            stream: SylowStream {
                stack: self.stream.stack.clone(),
                buffer: self.stream.buffer.clone(),
                tree: Arc::clone(&self.stream.tree),
            },
            splits: self.splits,
        }
    }
}

unsafe impl<S, const L: usize, C: SylowDecomposable<S> + Send, T: Send> Send for Seed<S, L, C, T> {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[derive(PartialEq, Eq, Debug)]
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
        SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_target(&[0, 1, 0])
            .into_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_leq_seq() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_flag(flags::LEQ)
            .add_target(&[2, 1, 0])
            .into_iter()
            .count();
        assert_eq!(count, 12);
    }

    #[test]
    pub fn test_generates_big_seq() {
        println!("A");
        let stream = SylowStreamBuilder::new()
            .add_target(&[0, 0, 0, 2, 0, 0, 0])
            .into_iter();
        let coords: Vec<SylowElem<Phantom, 7, FpNum<BIG_P>>> = stream.map(|(a, _)| a).collect();
        assert_eq!(coords.len(), 29 * 29 - 29);
        println!("B");

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>, ()>::new()
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
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new().add_target(&[0, 2, 1]);
        let stream_all = builder.into_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 2, 1])
            .into_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_seq() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_target(&[1, 0, 0])
            .add_target(&[0, 1, 0])
            .into_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
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
        let count = SylowStreamBuilder::<Phantom, 4, FpNum<13928643>, ()>::new()
            .add_flag(flags::LEQ)
            .add_target(&[0, 1, 1, 0])
            .into_iter()
            .count();

        assert_eq!(count, 91);
    }

    #[test]
    pub fn test_no_parabolic_seq() {
        let mut count = 0;
        SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_target(&[1, 0, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 1);
        let (mut x, _) = coords[0].clone();
        assert!(x != SylowElem::ONE);
        x = x.pow(2);
        assert!(x == SylowElem::ONE);

        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_target(&[0, 1, 0])
            .into_par_iter();
        assert_eq!(stream.count(), 2);
    }

    #[test]
    pub fn test_generates_big_par() {
        let stream = SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>, ()>::new()
            .add_target(&[0, 0, 0, 2, 0, 0, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 29 * 29 - 29);

        SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>, ()>::new()
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
        let builder = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new().add_target(&[0, 2, 1]);
        let stream_all = builder.into_par_iter();
        assert_eq!(stream_all.count(), 24);
    }

    #[test]
    pub fn test_skips_upper_half_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 2, 1])
            .into_par_iter();
        assert_eq!(stream.count(), 12);
    }

    #[test]
    pub fn test_multiple_targets_par() {
        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_target(&[1, 0, 0])
            .add_target(&[0, 1, 0])
            .into_par_iter();
        let coords = stream.collect::<Vec<_>>();
        assert_eq!(coords.len(), 3);

        let stream = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
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
        let coords = SylowStreamBuilder::<Phantom, 4, FpNum<13928643>, ()>::new()
            .add_flag(flags::LEQ)
            .add_target(&[0, 1, 1, 0])
            .into_par_iter()
            .collect::<Vec<_>>();

        assert_eq!(coords.len(), 91);
    }

    #[test]
    pub fn test_no_parabolic_par() {
        let count = AtomicUsize::new(0);
        SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
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
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_target(&[0, 1, 0])
            .add_target(&[0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 10);
    }

    #[test]
    pub fn test_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<61>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 4);
    }

    #[test]
    pub fn test_propagates_no_upper_half() {
        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[3, 1])
            .into_iter()
            .count();
        assert_eq!(count, 8);

        let count = SylowStreamBuilder::<Phantom, 2, FpNum<41>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[1, 1])
            .into_iter()
            .count();
        assert_eq!(count, 2);
    }

    #[test]
    pub fn test_quotient() {
        let res = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 3, 0])
            .set_quotient(Some([0, 2, 0]))
            .into_iter()
            .map(|(x, _)| x)
            .collect::<Vec<_>>();
        assert_eq!(res, vec![SylowElem::<Phantom, 3, FpNum<271>>::new([0, 1, 0])]);

        let count = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_target(&[0, 3, 0])
            .set_quotient(Some([0, 2, 0]))
            .into_iter()
            .count();
        assert_eq!(count, 2);

        let count = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_target(&[0, 3, 0])
            .set_quotient(Some([0, 1, 0]))
            .into_iter()
            .count();
        assert_eq!(count, 6);

        let count = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_PARABOLIC)
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[0, 3, 0])
            .set_quotient(Some([0, 1, 0]))
            .into_iter()
            .count();
        assert_eq!(count, 4);
    }

    #[test]
    pub fn test_generate_everything() {
        let count = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::LEQ)
            .add_target(&[1, 3, 1])
            .into_iter()
            .count();
        assert_eq!(count, 270);

        let count = SylowStreamBuilder::<Phantom, 3, FpNum<271>, ()>::new()
            .add_flag(flags::LEQ)
            .add_flag(flags::NO_UPPER_HALF)
            .add_target(&[1, 3, 1])
            .into_iter()
            .count();
        assert_eq!(count, 136);
    }
}
