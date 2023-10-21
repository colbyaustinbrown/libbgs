use crate::markoff::*;
use crate::numbers::*;
use crate::streams::*;

use rayon::iter::plumbing::*;
use rayon::iter::*;

/// A stream which can be run either in sequence or parallel, yielding Markoff numbers modulo `P`.
#[derive(Clone)]
pub struct CoordStream<'a, S, const L_HYPER: usize, const L_ELLIP: usize, const P: u128>
where
    FpNum<P>: SylowDecomposable<S>,
    QuadNum<P>: SylowDecomposable<S>,
{
    hyper_stream: Option<SylowSeqStream<S, L_HYPER, FpNum<P>>>,
    ellip_stream: Option<SylowSeqStream<S, L_ELLIP, QuadNum<P>>>,
    hyper_decomp: &'a SylowDecomp<S, L_HYPER, FpNum<P>>,
    ellip_decomp: &'a SylowDecomp<S, L_ELLIP, QuadNum<P>>,
}

impl<'a, S, const L_HYPER: usize, const L_ELLIP: usize, const P: u128>
    CoordStream<'a, S, L_HYPER, L_ELLIP, P>
where
    FpNum<P>: SylowDecomposable<S>,
    QuadNum<P>: SylowDecomposable<S>,
{
    /// Creates a new `CoordStream` with orders up to `limit`.
    pub fn new(
        hyper_decomp: &'a SylowDecomp<S, L_HYPER, FpNum<P>>,
        ellip_decomp: &'a SylowDecomp<S, L_ELLIP, QuadNum<P>>,
        limit: u128,
    ) -> CoordStream<'a, S, L_HYPER, L_ELLIP, P> {
        let hyper_stream = FactorStream::new(FpNum::FACTORS.prime_powers(), limit, true)
            .into_iter()
            .fold(
                SylowStreamBuilder::<S, L_HYPER, FpNum<P>>::new()
                    .add_flag(flags::NO_PARABOLIC)
                    .add_flag(flags::NO_UPPER_HALF)
                    .add_flag(flags::LEQ),
                |b, x| b.add_target(x),
            )
            .into_iter();
        let ellip_stream = FactorStream::new(QuadNum::FACTORS.prime_powers(), limit, true)
            .into_iter()
            .fold(
                SylowStreamBuilder::<S, L_ELLIP, QuadNum<P>>::new()
                    .add_flag(flags::NO_PARABOLIC)
                    .add_flag(flags::NO_UPPER_HALF)
                    .add_flag(flags::LEQ),
                |b, x| b.add_target(x),
            )
            .into_iter();
        CoordStream {
            hyper_stream: Some(hyper_stream),
            ellip_stream: Some(ellip_stream),
            hyper_decomp,
            ellip_decomp,
        }
    }

    /// Returns an iterator yielding pairs of coordinates without repeats up to permutation.
    pub fn upper_triangle(self) -> impl ParallelIterator<Item = (Coord<P>, Coord<P>)> + 'a
    where
        S: Clone + Send + Sync,
    {
        #[derive(Clone)]
        struct Helper<I>(I);
        impl<I: Iterator + Clone> Iterator for Helper<I> {
            type Item = (I, I::Item);
            fn next(&mut self) -> Option<Self::Item> {
                let iter = self.0.clone();
                let Some(a) = self.0.next() else {
                    return None;
                };
                Some((iter, a))
            }
        }
        Helper(self)
            .par_bridge()
            .flat_map(|(iter, a)| ParallelIterator::map(iter, move |b| (a.clone(), b)))
    }
}

impl<'a, S, const L_HYPER: usize, const L_ELLIP: usize, const P: u128> Iterator
    for CoordStream<'a, S, L_HYPER, L_ELLIP, P>
where
    FpNum<P>: SylowDecomposable<S>,
    QuadNum<P>: SylowDecomposable<S>,
{
    type Item = Coord<P>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(stream) = self.hyper_stream.as_mut() {
            if let Some(a) = stream.next() {
                return Some(Coord::from_chi_fp(&a, self.hyper_decomp));
            }
            self.hyper_stream = None;
        }
        if let Some(stream) = self.ellip_stream.as_mut() {
            if let Some(a) = stream.next() {
                return Some(Coord::from_chi_quad(&a, self.ellip_decomp));
            }
            self.ellip_stream = None;
        }
        None
    }
}

impl<'a, S, const L_HYPER: usize, const L_ELLIP: usize, const P: u128> ParallelIterator
    for CoordStream<'a, S, L_HYPER, L_ELLIP, P>
where
    S: Send + Sync,
    FpNum<P>: SylowDecomposable<S>,
    QuadNum<P>: SylowDecomposable<S>,
{
    type Item = Coord<P>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        let left = self.hyper_stream.map(|stream| {
            stream
                .parallelize()
                .map(|x| Coord::from_chi_fp(&x, self.hyper_decomp))
                .drive_unindexed(consumer.split_off_left())
        });
        let right = self.ellip_stream.map(|stream| {
            stream
                .parallelize()
                .map(|x| Coord::from_chi_quad(&x, self.ellip_decomp))
                .drive_unindexed(consumer.split_off_left())
        });
        match (left, right) {
            (Some(l), Some(r)) => consumer.to_reducer().reduce(l, r),
            (None, Some(r)) => r,
            (None, None) => consumer.into_folder().complete(),
            (Some(l), None) => l,
        }
    }
}
