use criterion::{criterion_group, criterion_main, Criterion};

use std::collections::HashSet;
use either::Either::*;

use libbgs::markoff::coord::Coord;
use libbgs::orbit_tester::*;
use libbgs::numbers::sylow_stream::*;
use libbgs::numbers::group::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::fp::FpStar;
use libbgs::numbers::quad_field::QuadField;
use libbgs::numbers::factorization::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

fn run_tester<'a>(stream: impl Iterator<Item = &'a u128>) {
    let mut tester = OrbitTester::new(&FpStar::<BIG_P> {});
    let mut count = 0;
    for x in stream {
        count += 1;
        tester = tester.add_target(*x);
    }

    tester.run();
}

fn criterion_benchmark(c: &mut Criterion) {
    let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
    let fp2_fact = Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]);
    let fp2 = QuadField::<BIG_P>::make();
    let fp_decomp = SylowDecomp::new(&FpStar::<BIG_P> {}, fp.clone());
    let fp2_decomp = SylowDecomp::new(&fp2, fp2_fact.clone());

    const LIMIT: u128 = 50;

    let mut fp_stream_builder = SylowStreamBuilder::new(&fp_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    let mut fp2_stream_builder = SylowStreamBuilder::new(&fp2_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);

    for d in fp.maximal_divisors(LIMIT) {
        fp_stream_builder = fp_stream_builder.add_target(d);
    }

    for d in fp2_fact.maximal_divisors(LIMIT) {
        fp2_stream_builder = fp2_stream_builder.add_target(d);
    }

    let stream = fp_stream_builder.build()
        .map(|x| Right((x.clone(), x)))
        .chain(fp2_stream_builder.build()
            .map(|x| Left((x.clone(), x)))
        )
        .map(|e| e.either(
            |(l, mut l_i)| { 
                l_i.invert(&fp2_decomp);
                Left((
                    l.to_product(&fp2_decomp), 
                    l_i.to_product(&fp2_decomp)
                ))
            },
            |(r, mut r_i)| { 
                r_i.invert(&fp_decomp);
                Right((
                    r.to_product(&fp_decomp), 
                    r_i.to_product(&fp_decomp)
                ))
            }
        ))
        .map(Coord::from_chi);

    let mut targets = HashSet::new();
    for x in stream {
        targets.insert(x.v());
    }

    c.bench_function("run orbit tester", |b| b.iter(|| run_tester(targets.iter())));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
