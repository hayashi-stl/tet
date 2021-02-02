use criterion::{black_box, criterion_group, criterion_main, Criterion, BatchSize};
use tet::Tets;

use nalgebra::Point3;
use rand::distributions::{Distribution, Uniform};
use rand::seq::SliceRandom;
use rand_distr::UnitSphere;
use rand_pcg::Pcg64;

type Pt3 = Point3<f64>;

const PCG_STATE: u128 = 0xcafef00dd15ea5e5;
const PCG_STREAM: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;

fn delaunay_tets_random(c: &mut Criterion) {
    let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
    let dist = Uniform::new_inclusive(-10.0, 10.0);
    let data = (0..10000)
        .map(|_| {
            let vals = dist.sample_iter(&mut rng).take(3).collect::<Vec<_>>();
            (Pt3::new(vals[0], vals[1], vals[2]), ())
        }).collect::<Vec<_>>();

    c.bench_function("delaunay_tets_random", |b| {
        b.iter_with_large_drop(|| Tets::<(), ()>::delaunay_from_vertices(data.iter().copied(), || ()))
    });
}

fn delaunay_tets_cospherical(c: &mut Criterion) {
    let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
    let data = UnitSphere.sample_iter(&mut rng).take(10000)
        .map(|vals| {
            (Pt3::new(vals[0], vals[1], vals[2]), ())
        }).collect::<Vec<_>>();

    c.bench_function("delaunay_tets_cospherical", |b| {
        b.iter_with_large_drop(|| Tets::<(), ()>::delaunay_from_vertices(data.iter().copied(), || ()))
    });
}

criterion_group!(bench_delaunay_tets, delaunay_tets_random, delaunay_tets_cospherical);
criterion_main!(bench_delaunay_tets);