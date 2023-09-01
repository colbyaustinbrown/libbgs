use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::*;

const P: u128 = 5109751;

#[derive(PartialEq, Eq)]
struct Phantom {}

impl Factored<Phantom, 4> for FpNum<P> {
    const FACTORS: Factorization<4> = Factorization::new([
        (2, 1),
        (3, 3),
        (5, 3),
        (757, 1)
    ]);
}

fn run_stream() {
    let g = SylowDecomp::<Phantom, 4, FpNum<P>>::new();

    let stream = SylowStreamBuilder::new(&g)
        .add_target(vec![0, 3, 2, 1])
        .build();

    let mut count = 0;
    for _ in stream {
        count += 1;
    }
    assert_eq!(count, 272160);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("run stream", |b| b.iter(|| run_stream()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
