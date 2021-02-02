use crate::{Pt3, VertexId, id_map::IdType};
use hilbert::{point_list, Point};
use float_ord::FloatOrd;

/// Sorts some points along the Hilbert curve.
/// This makes adjacent points in the list be close together in 3D.
pub(crate) fn hilbert_sort<V: Clone, I: IntoIterator<Item = (Pt3, V)>>(points: I) -> impl Iterator<Item = (VertexId, Pt3, V)> {
    let points = points.into_iter().collect::<Vec<_>>();

    // Range
    let min = points.iter().flat_map(|p| vec![p.0.x, p.0.y, p.0.z]).min_by_key(|x| FloatOrd(*x)).unwrap();
    let max = points.iter().flat_map(|p| vec![p.0.x, p.0.y, p.0.z]).max_by_key(|x| FloatOrd(*x)).unwrap();
    let max = max + max * f64::EPSILON;

    // Calculate scale to match Hilbert curve precision
    const PRECISION: usize = 32;

    let (mut h_points, bits) = point_list::make_points_f64(
        &points.iter().map(|p| [p.0.x, p.0.y, p.0.z]).collect::<Vec<_>>(),
        0, None, Some(PRECISION), 2.0f64.powi(PRECISION as i32) / (max - min)
    );

    // The good stuff: sorting
    Point::hilbert_sort(&mut h_points, bits);
    h_points.into_iter().map(move |p| { (VertexId(p.get_id() as IdType), points[p.get_id()].0, points[p.get_id()].1.clone()) })
}