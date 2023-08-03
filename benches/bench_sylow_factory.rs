use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::sylow_factory::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::factorization::*;

fn run_factory() {
        let fp = Factorization::new(vec![(2, 1), (3, 3), (5, 3), (757, 1)]);
        let g = SylowDecomp::new(&fp);

        let factory = SylowFactory::new(&g, vec![0, 3, 2, 1], flags::NONE); 

        let mut count = 0;
        for y in factory { count += 1; }
        assert_eq!(count, 272160);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("run factory", |b| b.iter(|| run_factory()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

