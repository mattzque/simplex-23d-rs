use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simplex_23d::Simplex;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("sample2d", |b| {
        b.iter_with_setup(
            || Simplex::new(42),
            |noise| noise.sample2d(black_box(1.0), black_box(1.0)),
        )
    });
    c.bench_function("sample3d", |b| {
        b.iter_with_setup(
            || Simplex::new(42),
            |noise| noise.sample3d(black_box(1.0), black_box(1.0), black_box(1.0)),
        )
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
