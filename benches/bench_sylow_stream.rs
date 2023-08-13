use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::fp::*;
use libbgs::numbers::sylow_stream::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::factorization::*;

fn run_stream() {
        let fp = FpStar::<5109751> {};
        let g = SylowDecomp::new(&fp, Factorization::new(vec![(2, 1), (3, 3), (5, 3), (757, 1)]));

        let stream = SylowStreamBuilder::new(&g)
            .add_target(vec![0, 3, 2, 1]) 
            .build();

        let mut count = 0;
        for _ in stream { count += 1; }
        assert_eq!(count, 272160);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("run stream", |b| b.iter(|| run_stream()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

