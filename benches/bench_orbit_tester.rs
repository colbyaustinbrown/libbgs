use criterion::{criterion_group, criterion_main, Criterion};

use std::collections::HashSet;

use rayon::iter::*;

use libbgs::markoff::*;
use libbgs::numbers::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

#[derive(PartialEq, Eq)]
struct Phantom {}

impl Factor<Phantom, 7> for FpNum<BIG_P> {
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

impl Factor<Phantom, 11> for QuadNum<BIG_P> {
    const FACTORS: Factorization<11> = Factorization::new([
        (2, 4),
        (3, 1),
        (5, 2),
        (11, 2),
        (17, 1),
        (19, 1),
        (23, 1),
        (97, 1),
        (757, 1),
        (1453, 1),
        (8689, 1),
    ]);
}

fn run_tester<'a>(stream: impl Iterator<Item = &'a u128>) {
    let mut tester = OrbitTester::<BIG_P>::new();
    for x in stream {
        tester.add_target(*x);
    }

    tester.run();
}

fn criterion_benchmark(c: &mut Criterion) {
    const LIMIT: u128 = 50;
    let fp_decomp = SylowDecomp::new();
    let fp2_decomp = SylowDecomp::new();

    let mut fp_stream_builder = SylowStreamBuilder::new()
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    let mut fp2_stream_builder = SylowStreamBuilder::new()
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);

    for d in FpNum::FACTORS.maximal_divisors(LIMIT) {
        fp_stream_builder = fp_stream_builder.add_target(d);
    }

    for d in QuadNum::FACTORS.maximal_divisors(LIMIT) {
        fp2_stream_builder = fp2_stream_builder.add_target(d);
    }

    let targets: HashSet<_> = fp_stream_builder
        .into_par_iter()
        .map(|x| Coord::from_chi_fp(&x, &fp_decomp))
        .chain(
            fp2_stream_builder
                .into_par_iter()
                .map(|x| Coord::from_chi_quad(&x, &fp2_decomp)),
        )
        .map(|x| x.into())
        .collect();

    let mut group = c.benchmark_group("orbits");
    group.sample_size(20);
    group.bench_function("run orbit tester", |b| {
        b.iter(|| run_tester(targets.iter()))
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
