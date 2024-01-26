#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use rayon::iter::*;

use libbgs::markoff::*;
use libbgs::numbers::*;
use libbgs::streams::*;
use libbgs::util::*;

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
    coset_max: AtomicUsize,
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
            RotOrder::Parabola => false,
            RotOrder::Hyperbola(ord) => ord <= self.hyper_lim,
            RotOrder::Ellipse(ord) => ord <= self.ellip_lim,
        }
    }
}

#[derive(Debug)]
enum Check<const L: usize> {
    Cosets([usize; L]),
    SmallOrders(u128),
    // Skip,
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
        coset_max: AtomicUsize::new(0),
    };

    // Magic number used to permute cosets of <chi> to ensure all (s*chi + (s*chi)^-1) have order
    // dividing 2(p - 1) and not dividing (p - 1)
    let magic = (1..P*P)
        .map(|i| {
            let j = standard_affine_shift(P * 2, i);
            QuadNum::<P>::steinitz(j)
        })
        .filter(|c| *c != QuadNum::ZERO)
        .find_map(|c| {
            let twos = (QuadNum::<P>::FACTORS[0].1 + 1) as u128;
            let pow = (P*P - 1) / intpow::<0>(2, twos);
            let res = c.pow(pow);
            if res.pow(intpow::<0>(2, twos - 1)) == QuadNum::ONE {
                None 
            } else {
                Some(res)
            }
        })
        .unwrap();
    let magic = if QuadNum::<P>::FACTORS[0].1 == 1 {
        magic * QuadNum::find_sylow_generator(1)
    } else {
        magic
    };

    let (a, b) = rayon::join(
        || process_trie(&elements_count, hyper_lim, &ctx.hyper_decomp, &ctx, |k, s| { k * (s + s.inverse())}),
        || process_trie(&elements_count, ellip_lim, &ctx.ellip_decomp, &ctx, |k, s| { 
                let fix = s * magic;
                let b = fix + fix.inverse();
                assert_eq!(b.0, FpNum::ZERO);
                k * b.1
            }),
    );
    let dur = now.elapsed();
    println!(
        "{} {} {} {} {:?} {} {} {}",
        P,
        dur.as_millis(),
        hyper_endgame,
        ellip_endgame,
        middle_game,
        ctx.coset_max.into_inner(),
        a,
        b
    );
}

fn process_trie<const P: u128, C>(
    counts: &HashMap<u128, u128>,
    limit: u128,
    decomp: &SylowDecomp<Ph, { C::LENGTH }, C>,
    ctx: &Context<P>,
    get_coset_repr: impl Fn(FpNum<P>, C) -> FpNum<P> + Send + Sync,
) -> u128
where
    C: SylowDecomposable<Ph> + FromChi<Ph, P> + Send + Sync + Copy + std::fmt::Debug,
    FpNum<P>: Factor<Ph>,
    QuadNum<P>: Factor<Ph>,
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
            None => Check::Cosets(*ds),
        }
    });
    SylowStreamBuilder::new_with_trie(&trie)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::LEQ)
        .add_targets_leq(limit)
        .into_par_iter()
        .map(|(chi, check)| {
            let a = C::from_chi(&chi, &decomp);
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
                    let chi_conj = C::from_chi_conj(&chi, &decomp).inverse();

                    SylowStreamBuilder::new_with_trie(&trie)
                    .add_flag(flags::NO_UPPER_HALF)
                    .add_targets_leq(P + 1)
                    .set_quotient(Some(*gen))
                    .into_par_iter()
                    .map(|(x, _)| {
                        let b = a * get_coset_repr(chi_conj, x.to_product(decomp));

                        if a == FpNum::from(0) && b == FpNum::from(0) {
                            return 0;
                        }
                        if !ctx.is_small(&Coord(b)) {
                            return 0;
                        }
                        let Some(mut it) = Coord(a).part(Coord(b)) else {
                            panic!("Attempted to look at coset solutions that don't exist: P={} a={} b={}.", P, u128::from(a), u128::from(b));
                        };
                        let mut count = 0;
                        let res = if it
                            .take(50)
                            .all(|c| {
                                count += 1;
                                ctx.is_small(&c)
                            })
                        {
                            chi.order()
                        } else {
                            0
                        };
                        ctx.coset_max.fetch_max(count, Ordering::Relaxed);
                        res
                    })
                    .sum()
                },
                // Check::Skip => 0,
            }
        })
        .sum()
}

impl_factors!(Ph, 4000..5000);

macro_rules! go {
    ($($P:literal),+$(,)?) => {$(process::<$P>();)+};
}

fn main() {
    primes!(go, 4000..5000);
}
