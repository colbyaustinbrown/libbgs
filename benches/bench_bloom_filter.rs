use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::markoff::*;

fn query(filter: &BloomFilter<u64>) {
    for i in 1..10_000 {
        filter.is_member_prob(&i); 
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let hashes: Vec<Box<dyn Fn(&u64) -> usize>> = vec![
        Box::new(|x| (x % 10_000) as usize),
        Box::new(|x| (((x % 997) * (x % 991)) % 10_000) as usize),
    ];
    let mut filter = BloomFilter::<u64>::new(10_000, hashes);
    for i in 50_000..51_000 {
        filter.add(&i);
    }
    
    c.bench_function("query bloom filter", |b| b.iter(|| query(&filter)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
