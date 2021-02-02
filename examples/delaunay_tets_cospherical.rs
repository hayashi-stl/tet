//! Tetrahedralize a random point set

use rand_distr::UnitSphere;
use tet::Tets;

use nalgebra::Point3;
use rand::distributions::{Distribution, Uniform};
use rand_pcg::Pcg64;

type Pt3 = Point3<f64>;

const PCG_STATE: u128 = 0xcafef00dd15ea5e5;
const PCG_STREAM: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;

fn main() {
    let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
    let data = UnitSphere.sample_iter(&mut rng).take(10000)
        .map(|vals| {
            (Pt3::new(vals[0], vals[1], vals[2]), ())
        }).collect::<Vec<_>>();

    let mesh = Tets::<(), ()>::delaunay_from_vertices(data, || ());
    println!("Num tets including ghosts: {:?}", mesh.num_tets());
}