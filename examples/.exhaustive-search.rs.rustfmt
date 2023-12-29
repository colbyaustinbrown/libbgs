// Requires nightly!
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use rayon::iter::*;

use libbgs::markoff::*;
use libbgs::numbers::*;
use libbgs::streams::*;

#[derive(Debug)]
struct Ph {}

struct Context<const P: u128>
where
    FpNum<P>: Factor<Ph>,
    QuadNum<P>: Factor<Ph>,
    [(); FpNum::<P>::LENGTH]:,
    [(); QuadNum::<P>::LENGTH]:,
{
    hyper_lim: u128,
    ellip_lim: u128,
    hyper_decomp: SylowDecomp<Ph, { FpNum::<P>::LENGTH }, FpNum<P>>,
    ellip_decomp: SylowDecomp<Ph, { QuadNum::<P>::LENGTH }, QuadNum<P>>,
    coset_searches: AtomicUsize,
}

impl<const P: u128> Context<P>
where
    FpNum<P>: Factor<Ph>,
    QuadNum<P>: Factor<Ph>,
    [(); FpNum::<P>::LENGTH]:,
    [(); QuadNum::<P>::LENGTH]:,
{
    fn is_small(&self, c: &Coord<P>) -> bool {
        match c.rot_order() {
            (_, Conic::Parabola) => false,
            (ord, Conic::Hyperbola) => ord < self.hyper_lim,
            (ord, Conic::Ellipse) => ord < self.ellip_lim,
        }
    }
}

#[derive(Debug)]
enum Check<const L: usize> {
    Cosets([usize; L]),
    SmallOrders(u128),
    Skip,
}

fn process<const P: u128>()
where
    FpNum<P>: Factor<Ph>,
    QuadNum<P>: Factor<Ph>,
    [(); FpNum::<P>::LENGTH]:,
    [(); QuadNum::<P>::LENGTH]:,
{
    let now = Instant::now();

    let hyper_decomp = SylowDecomp::<Ph, { FpNum::<P>::LENGTH }, FpNum<P>>::new();
    let ellip_decomp = SylowDecomp::<Ph, { QuadNum::<P>::LENGTH }, QuadNum<P>>::new();

    let (hyper_endgame, ellip_endgame) = Coord::<P>::endgame();

    let mut factors = DivisorStream::new(FpNum::<P>::FACTORS.factors(), hyper_endgame, false)
        .map(|x| FpNum::<P>::FACTORS.from_powers(&x))
        .chain(
            DivisorStream::new(QuadNum::<P>::FACTORS.factors(), ellip_endgame, false)
                .map(|x| QuadNum::<P>::FACTORS.from_powers(&x)),
        )
        .collect::<Vec<_>>();
    factors.sort_unstable();

    let mut elements_count = HashMap::<u128, u128>::new();
    let mut middle_game = None;

    for t in factors {
        let mut rhs = 0.0;
        let mut count = 0;
        for (d, c) in DivisorStream::new(FpNum::<P>::FACTORS.factors(), t, true)
            .map(|x| {
                (
                    FpNum::<P>::FACTORS.from_powers(&x),
                    FpNum::<P>::count_elements_of_order(&x),
                )
            })
            .chain(
                DivisorStream::new(QuadNum::<P>::FACTORS.factors(), t, true).map(|x| {
                    (
                        QuadNum::<P>::FACTORS.from_powers(&x),
                        QuadNum::<P>::count_elements_of_order(&x),
                    )
                }),
            )
        {
            let a = (6.0 * (t as f64) * (d as f64)).powf(1.0 / 3.0);
            let b = 4.0 * (t as f64) * (d as f64) / (P as f64);
            rhs += 1.5 * (if a >= b { a } else { b });
            count += c;
        }
        if (t as f64) >= rhs {
            middle_game.get_or_insert(t);
        } else {
            middle_game = None;
        }
        elements_count.insert(t, count);
    }
    let (hyper_lim, ellip_lim) = if let Some(min) = middle_game {
        (
            std::cmp::min(min, hyper_endgame),
            std::cmp::min(min, ellip_endgame),
        )
    } else {
        (hyper_endgame, ellip_endgame)
    };

    let ctx = Context {
        hyper_decomp,
        ellip_decomp,
        hyper_lim,
        ellip_lim,
        coset_searches: AtomicUsize::new(0),
    };

    let (a, b) = rayon::join(
        || process_trie(&elements_count, hyper_lim, &ctx.hyper_decomp, &ctx),
        || process_trie(&elements_count, ellip_lim, &ctx.ellip_decomp, &ctx),
    );
    let dur = now.elapsed();
    println!(
        "{} {} {} {} {:?} {} {} {}",
        P,
        dur.as_millis(),
        hyper_endgame,
        ellip_endgame,
        middle_game,
        ctx.coset_searches.into_inner(),
        a,
        b
    );
}

fn process_trie<const P: u128, C>(
    counts: &HashMap<u128, u128>,
    limit: u128,
    decomp: &SylowDecomp<Ph, { C::LENGTH }, C>,
    ctx: &Context<P>,
) -> u128
where
    C: SylowDecomposable<Ph> + Send + Sync,
    FpNum<P>: Factor<Ph>,
    QuadNum<P>: Factor<Ph>,
    Coord<P>: FromChi<Ph, P, C>,
    [(); C::LENGTH]:,
    [(); FpNum::<P>::LENGTH]:,
    [(); QuadNum::<P>::LENGTH]:,
{
    let trie = FactorTrie::<Ph, { C::LENGTH }, C, Check<{ C::LENGTH }>>::new_with(|ds, _| {
        let ord = C::FACTORS.from_powers(ds);
        let cosets = C::SIZE / ord;
        match counts.get(&ord) {
            Some(count) if limit == C::SIZE - 1 || *count > cosets => Check::Cosets(*ds),
            Some(_) => Check::SmallOrders(ord),
            None => Check::Skip,
        }
    });
    SylowStreamBuilder::new_with_trie(&trie)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::NO_UPPER_HALF)
        .add_targets_leq(limit)
        .into_par_iter()
        .map(|(chi, check)| {
            let a = Coord::from_chi(&chi, &decomp).0;
            match check {
                Check::SmallOrders(ord) => {
                    let it = CoordStream::new(&ctx.hyper_decomp, &ctx.ellip_decomp, *ord, *ord)
                        .into_par_iter();
                    ParallelIterator::flat_map(it, |b| {
                        let b = b.0;
                        match (a * a * b * b - 4 * (a * a + b * b)).int_sqrt() {
                            None => Vec::new(),
                            Some(disc) if disc == FpNum::<P>::ZERO => {
                                vec![a * b * FpNum::TWO_INV]
                            },
                            Some(disc) => {
                                vec![(a * b - disc) * FpNum::TWO_INV, (a * b + disc) * FpNum::TWO_INV]
                            },
                        }
                    })
                    .filter(|c| ctx.is_small(&Coord(*c)))
                    .count() as u128
                }
                Check::Cosets(gen) => {
                    let k = a * Coord::from_chi_conj(&chi, &decomp).0.inverse();
                    SylowStreamBuilder::new_with_trie(&trie)
                    .add_flag(flags::NO_UPPER_HALF)
                    .add_targets_leq(P + 1)
                    .set_quotient(Some(*gen))
                    .into_par_iter()
                    .map(|x| {
                        ctx.coset_searches.fetch_add(1, Ordering::Relaxed);
                        let b = Coord(k * Coord::from_chi(&x.0, &decomp).0);
                        if b.rot_order().0 < limit {
                            return 0;
                        }
                        let Some(mut it) = Coord(a).part(b) else {
                            panic!("Attempted to look at coset solutions that don't exist: {} {} {}.", P, u128::from(a), u128::from(b));
                        };
                        if it.all(|c| ctx.is_small(&c)) {
                            x.0.order()
                        } else {
                            0
                        }
                    })
                    .sum()
                },
                Check::Skip => 0,
            }
        })
        .sum()
}

impl_factors!(Ph, 1000..2000);

macro_rules! go {
    ($($P:literal),+$(,)?) => {$(process::<$P>();)+};
}

fn main() {
    primes!(go, 1000..2000);
}
