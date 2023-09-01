use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

fn find_sqrt() {
    for i in 3..103 {
        let mut x = QuadNum::<BIG_P>::int_sqrt(i);
        let y = x.clone();
        assert_ne!(x, i);
        x = x.multiply(&y);
        assert_eq!(x, i);
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("find sqrt", |b| b.iter(|| find_sqrt()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
