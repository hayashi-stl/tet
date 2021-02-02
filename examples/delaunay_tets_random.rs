//! Tetrahedralize a random point set

use tet::Tets;

use nalgebra::Point3;
use rand::distributions::{Distribution, Uniform};
use rand_pcg::Pcg64;

type Pt3 = Point3<f64>;

const PCG_STATE: u128 = 0xcafef00dd15ea5e5;
const PCG_STREAM: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;

fn main() {
    let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
    let dist = Uniform::new_inclusive(-10.0, 10.0);
    let data = (0..10000)
        .map(|_| {
            let vals = dist.sample_iter(&mut rng).take(3).collect::<Vec<_>>();
            (Pt3::new(vals[0], vals[1], vals[2]), ())
        });

    let mesh = Tets::<(), ()>::delaunay_from_vertices(data, || ());
    println!("Num tets including ghosts: {:?}", mesh.num_tets());
}