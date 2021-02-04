//! Intersection testing for PLC elements

use simplicity as sim;
use std::iter::FromIterator;

use nalgebra::{Isometry3, Unit};
use ncollide3d::bounding_volume::{self, AABB};
use ncollide3d::query::{self, Ray};
use ncollide3d::shape::{Plane, Segment};
use ncollide3d::{
    bounding_volume::BoundingVolume,
    partitioning::{SimultaneousVisitor, VisitStatus},
    query::{PointQuery, RayCast},
    shape::Shape,
};

use crate::{Pt3, Vec1, Vec3};

use super::{EdgeId, Element, FaceId, Plc, Vertex, VertexId};

fn point_in_face<V, E, F>(
    plc: &Plc<V, E, F>,
    point: Pt3,
    face: Vec<Vec<(VertexId, &Vertex<V>)>>,
    center: Pt3,
    normal: Unit<Vec3>,
) -> bool {
    // Transform face to 2D
    let mut tangent = normal.cross(&Vec3::new(1.0, 0.0, 0.0));
    if normal.x.abs() > 0.9 {
        tangent = normal.cross(&Vec3::new(0.0, 1.0, 0.0));
    }
    let tangent = Unit::new_normalize(tangent);

    let transform = Isometry3::look_at_lh(&center, &(center + normal.xyz()), &tangent);
    let point = transform * point;
    // Check each polygon
    face.into_iter()
        .map(|ring| {
            let mut ring = ring
                .into_iter()
                .map(|(_, v)| transform * v.position())
                .collect::<Vec<_>>();
            ring.push(ring[0]);
            // Calculate intersection count
            ring.windows(2)
                .map(|vertices| {
                    // Technique taken from Simulation of Simplicity.
                    let mut pts = [point.xy(), vertices[0].xy(), vertices[1].xy()];
                    if pts[1].y > pts[2].y {
                        pts.swap(1, 2);
                    }

                    (pts[1].y <= pts[0].y
                        && pts[0].y < pts[2].y
                        && sim::orient_2d(&pts, |l, i| l[i].coords, 0, 1, 2))
                        as usize
                })
                .sum::<usize>()
                % 2
        })
        .sum::<usize>()
        % 2
        != 0
}

pub(super) struct IntersectionTest<'a, V, E, F> {
    pub(super) plc: &'a Plc<V, E, F>,
    pub(super) tolerance: f64,
    pub(super) intersections: Vec<(Element, Element, f64)>,
}

impl<'a, V, E, F> SimultaneousVisitor<Element, AABB<f64>> for IntersectionTest<'a, V, E, F> {
    fn visit(
        &mut self,
        left_bv: &AABB<f64>,
        left_data: Option<&Element>,
        right_bv: &AABB<f64>,
        right_data: Option<&Element>,
    ) -> VisitStatus {
        if !left_bv.intersects(right_bv) {
            return VisitStatus::Stop;
        }

        let (left, right) = match (left_data, right_data) {
            (Some(left), Some(right)) => (left, right),
            _ => return VisitStatus::Continue,
        };

        // Actual intersection
        // Do NOT check for intersections if the elements share a vertex,
        // as that means that either a lower-order intersection would catch the problem anyway
        // or the elements are allowed to intersect because the intersection is in the PLC.
        match (left, right) {
            (Element::Vertex(v0), Element::Vertex(v1)) => {
                if v0 < v1 {
                    let dist_squared =
                        (self.plc[*v0].position() - self.plc[*v1].position()).norm_squared();
                    if dist_squared < self.tolerance * self.tolerance {
                        self.intersections.push((
                            Element::Vertex(*v0),
                            Element::Vertex(*v1),
                            dist_squared.sqrt(),
                        ));
                    }
                }
            }

            (Element::Vertex(v), Element::Edge(e)) => {
                let ev = self.plc[*e].vertices();
                if !ev.contains(v) {
                    let edge = Segment::new(self.plc[ev[0]].position(), self.plc[ev[1]].position());
                    let dist = edge.distance_to_point(
                        &Isometry3::identity(),
                        &self.plc[*v].position(),
                        true,
                    );
                    if dist < self.tolerance {
                        self.intersections
                            .push((Element::Vertex(*v), Element::Edge(*e), dist));
                    }
                }
            }

            (Element::Edge(e0), Element::Edge(e1)) => {
                let ev0 = self.plc[*e0].vertices();
                let ev1 = self.plc[*e1].vertices();
                if ev0 < ev1 && !ev0.contains(&ev1[0]) && !ev0.contains(&ev1[1]) {
                    let s0 = Segment::new(self.plc[ev0[0]].position(), self.plc[ev0[1]].position());
                    let s1 = Segment::new(self.plc[ev1[0]].position(), self.plc[ev1[1]].position());
                    let dist =
                        query::distance(&Isometry3::identity(), &s0, &Isometry3::identity(), &s1);
                    if dist < self.tolerance {
                        self.intersections
                            .push((Element::Edge(*e0), Element::Edge(*e1), dist));
                    }
                }
            }

            (Element::Vertex(v), Element::Face(f)) => {
                let fv = self.plc[*f]
                    .vertices(self.plc)
                    .map(Vec::from_iter)
                    .collect::<Vec<_>>();
                if !fv.iter().flatten().any(|(vertex, _)| v == vertex) {
                    let plane = Plane::new(self.plc[*f].plane.1);
                    let center = self.plc[*f].plane.0;
                    let p0 = self.plc[*v].position();

                    // Why does Isometry3::translation take 3 separate arguments?!
                    // Also, this distance is negative if the point is "inside" the plane.
                    let dist = plane
                        .distance_to_point(
                            &Isometry3::translation(center.x, center.y, center.z),
                            &p0,
                            false,
                        )
                        .abs();
                    if dist < self.tolerance
                        && point_in_face(self.plc, p0, fv, center, plane.normal)
                    {
                        self.intersections
                            .push((Element::Vertex(*v), Element::Face(*f), dist));
                    }
                }
            }

            (Element::Edge(e), Element::Face(f)) => {
                let fv = self.plc[*f]
                    .vertices(self.plc)
                    .map(Vec::from_iter)
                    .collect::<Vec<_>>();
                let ev = self.plc[*e].vertices();
                if !fv
                    .iter()
                    .flatten()
                    .any(|(vertex, _)| ev[0] == *vertex || ev[1] == *vertex)
                {
                    let plane = Plane::new(self.plc[*f].plane.1);
                    let (center, normal) = self.plc[*f].plane;
                    let positions = [self.plc[ev[0]].position(), self.plc[ev[1]].position()];
                    let ray = Ray::new(positions[0], positions[1] - positions[0]);

                    // Why does Isometry3::translation take 3 separate arguments?!
                    if let Some(toi) = plane.toi_with_ray(
                        &Isometry3::translation(center.x, center.y, center.z),
                        &ray,
                        1.0,
                        false,
                    ) {
                        let point = positions[0] + (positions[1] - positions[0]) * toi;
                        if point_in_face(self.plc, point, fv, center, normal) {
                            self.intersections
                                .push((Element::Edge(*e), Element::Face(*f), 0.0));
                        }
                    }
                }
            }

            // Intersection/closeness implied by reverse direction
            (Element::Face(_), Element::Edge(_))
            | (Element::Face(_), Element::Vertex(_))
            | (Element::Edge(_), Element::Vertex(_)) => {}

            // Intersection implied by lower-order (face-edge) intersection
            (Element::Face(_), Element::Face(_)) => {}
        }

        VisitStatus::Continue
    }
}

pub(super) fn vertex_bound<V, E, F>(
    plc: &Plc<V, E, F>,
    vertex: VertexId,
    tolerance: f64,
) -> (Element, AABB<f64>) {
    let pos = plc[vertex].position();
    let aabb = AABB::new(
        pos - Vec1::new(tolerance).xxx(),
        pos + Vec1::new(tolerance).xxx(),
    );
    (Element::Vertex(vertex), aabb)
}

pub(super) fn edge_bound<V, E, F>(
    plc: &Plc<V, E, F>,
    edge: EdgeId,
    tolerance: f64,
) -> (Element, AABB<f64>) {
    let vs = plc[edge].vertices();
    let segment = Segment::new(plc[vs[0]].position(), plc[vs[1]].position());
    let mut aabb = segment.aabb(&Isometry3::identity());
    aabb.mins -= Vec1::new(tolerance).xxx();
    aabb.maxs += Vec1::new(tolerance).xxx();
    (Element::Edge(edge), aabb)
}

pub(super) fn face_bound<V, E, F>(
    plc: &Plc<V, E, F>,
    face: FaceId,
    tolerance: f64,
) -> (Element, AABB<f64>) {
    let pts = plc[face].vertices(plc).flatten().map(|(_, v)| &v.position);
    let mut aabb = bounding_volume::local_point_cloud_aabb(pts);
    aabb.mins -= Vec1::new(tolerance).xxx();
    aabb.maxs += Vec1::new(tolerance).xxx();
    (Element::Face(face), aabb)
}
