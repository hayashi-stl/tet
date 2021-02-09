use nalgebra::Isometry3;

use super::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeIntersectionRaw {
    pub other: EdgeOrTri,
    /// Unnormalized barycentric coordinates of intersection point on edge.
    bary: [f64; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeIntersection {
    pub other: EdgeOrTri,
    pub point: Pt3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriIntersectionRaw {
    pub edge: [VertexId; 2],
    /// Unnormalized barycentric coordinates of intersection point on triangle.
    bary: [f64; 3],
}

impl TriIntersection {
    fn new(edge: [VertexId; 2], point: Pt3) -> Self {
        Self { edge, point }
    }

    fn sorted(self) -> Self {
        Self { edge: sorted_2(self.edge), point: self.point }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriIntersection {
    pub edge: [VertexId; 2],
    pub point: Pt3,
}

impl EdgeIntersection {
    fn new(other: EdgeOrTri, point: Pt3) -> Self {
        Self { other, point }
    }

    fn sorted(self) -> Self {
        Self { other: self.other.sorted(), point: self.point }
    }
}

impl<V, T> TetMesh<V, T> {
    /// Gets the normal of a triangle. This might not be normalized.
    pub fn tri_normal(&self, v0: VertexId, v1: VertexId, v2: VertexId) -> Vec3 {
        let p0 = self[v0].position().coords;
        let p1 = self[v1].position().coords;
        let p2 = self[v2].position().coords;

        Vec3::new(
            robust_geo::orient_2d(p0.yz(), p1.yz(), p2.yz()),
            robust_geo::orient_2d(p0.zx(), p1.zx(), p2.zx()),
            robust_geo::orient_2d(p0.xy(), p1.xy(), p2.xy()),
        )
    }

    /// Gets the transform needed to put the triangle onto the xy plane.
    pub fn tri_xy(&self, v0: VertexId, v1: VertexId, v2: VertexId) -> Isometry3<f64> {
        let center = self[v0].position();
        let tangent = self[v1].position() - center;
        let normal = self.tri_normal(v0, v1, v2);
        Isometry3::look_at_lh(&center, &(center + normal), &tangent)
    }

    /// Gets whether these points are oriented positive.
    /// Any vertex is allowed to be the ghost vertex and this should just work™.
    pub fn orient_3d(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId) -> f64 {
        let p0 = self.idx_ori(v0);
        let p1 = self.idx_ori(v1);
        let p2 = self.idx_ori(v2);
        let p3 = self.idx_ori(v3);
        (if [v0, v1, v2, v3].contains(&Self::GHOST) {-1.0} else {1.0}) * robust_geo::orient_3d(p0, p1, p2, p3)
    }

    /// Gets whether these points are oriented positive.
    /// Uses simulation of simplicity to avoid ties.
    /// Any vertex is allowed to be the ghost vertex and this should just work™.
    pub fn orient_3d_sim(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId) -> bool {
        [v0, v1, v2, v3].contains(&Self::GHOST)
            != sim::orient_3d(self, Self::idx_ori, v0, v1, v2, v3)
    }

    /// Gets the quality of a tet, measured as smallest distance between two opposite simplexes.
    pub fn quality(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId) -> f64 {
        if [v0, v1, v2, v3].contains(&Self::GHOST) {
            if self.center.x.is_nan() || !sim::orient_3d(self, Self::idx_ori, v0, v1, v2, v3) {
                f64::INFINITY
            } else {
                -1.0
            }
        } else {
            // Make sure the number returned is the same given the same tet, no matter the (even) order of the vertices.
            let [v0, v1, v2, v3] = even_sorted_4([v0, v1, v2, v3]);

            let p0 = self.idx_ori(v0);
            let p1 = self.idx_ori(v1);
            let p2 = self.idx_ori(v2);
            let p3 = self.idx_ori(v3);
            let volx6 = robust_geo::orient_3d(p0, p1, p2, p3);
            if volx6 <= 0.0 {
                return -1.0;
            }

            //{
            //    let mut vs = [v0.0, v1.0, v2.0, v3.0];
            //    vs.sort();
            //    if vs == [1005, 1006, 1007, 1008] {
            //        println!("Vertices: {} {} {} {}", v0.0, v1.0, v2.0, v3.0);
            //        println!("Vol x6: {}", volx6);
            //        println!("Tri 0 cross: {}", (p2 - p1).cross(&(p3 - p1)).norm());
            //        println!("Tri 1 cross: {}", (p3 - p2).cross(&(p0 - p2)).norm());
            //        println!("Tri 2 cross: {}", (p0 - p3).cross(&(p1 - p3)).norm());
            //        println!("Tri 3 cross: {}", (p1 - p0).cross(&(p2 - p0)).norm());
            //        println!("Edge 01 cross: {}", (p1 - p0).cross(&(p3 - p2)).norm());
            //        println!("Edge 02 cross: {}", (p2 - p0).cross(&(p3 - p1)).norm());
            //        println!("Edge 03 cross: {}", (p3 - p0).cross(&(p1 - p2)).norm());
            //    }
            //}

            let max_area = [[p0, p1, p0, p2], [p3, p2, p3, p1], [p2, p3, p2, p0], [p1, p0, p1, p3],
                [p0, p1, p2, p3], [p0, p2, p1, p3], [p0, p3, p1, p2]].iter().map(|[pa, pb, pc, pd]|
                    (pb - pa).cross(&(pd - pc)).norm())
                    .max_by_key(|area| FloatOrd(*area)).unwrap();
            volx6 / max_area
        }
    }

    /// Checks whether these points are oriented positive and do not form a sliver.
    pub fn non_sliver(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId) -> bool {
        self.quality(v0, v1, v2, v3) >= self.tolerance
    }

    /// Gets whether the last point is in the circumsphere of the first 4 points.
    /// Returns 0 in case of a tie.
    /// [v0, v1, v2, v3] must have positive orientation.
    /// Assumes that none of the points equal and that the last point isn't the ghost vertex.
    pub fn in_sphere(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> f64 {
        let p0 = self.idx_ori(v0);
        let p1 = self.idx_ori(v1);
        let p2 = self.idx_ori(v2);
        let p3 = self.idx_ori(v3);
        let p4 = self[v4].position().coords;
        if v0 == Self::GHOST {
            robust_geo::orient_3d(p3, p2, p1, p4)
        } else if v1 == Self::GHOST {
            robust_geo::orient_3d(p2, p3, p0, p4)
        } else if v2 == Self::GHOST {
            robust_geo::orient_3d(p1, p0, p3, p4)
        } else if v3 == Self::GHOST {
            robust_geo::orient_3d(p0, p1, p2, p4)
        } else {
            robust_geo::in_sphere(p0, p1, p2, p3, p4)
        }
    }

    /// Gets whether the last point is in the circumsphere of the first 4 points.
    /// Uses simulation of simplicity to avoid ties.
    /// [v0, v1, v2, v3] must have positive orientation.
    /// Assumes that none of the points equal and that the last point isn't the ghost vertex.
    //pub fn in_sphere_sim(
    //    &self,
    //    v0: VertexId,
    //    v1: VertexId,
    //    v2: VertexId,
    //    v3: VertexId,
    //    v4: VertexId,
    //) -> bool {
    //    if v0 == Self::GHOST {
    //        sim::orient_3d(self, Self::idx, v3, v2, v1, v4)
    //    } else if v1 == Self::GHOST {
    //        sim::orient_3d(self, Self::idx, v2, v3, v0, v4)
    //    } else if v2 == Self::GHOST {
    //        sim::orient_3d(self, Self::idx, v1, v0, v3, v4)
    //    } else if v3 == Self::GHOST {
    //        sim::orient_3d(self, Self::idx, v0, v1, v2, v4)
    //    } else {
    //        sim::in_sphere(self, Self::idx, v0, v1, v2, v3, v4)
    //    }
    //}

    /// Gets whether the last point is close to the tet formed by the first 4 points using the tolerance.
    /// The last point can't be the ghost vertex.
    /// The tet must be oriented positive.
    pub fn tet_close_to_vertex(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
        let p0 = self.idx_ori(v0);
        let p1 = self.idx_ori(v1);
        let p2 = self.idx_ori(v2);
        let p3 = self.idx_ori(v3);
        let p4 = self.idx_ori(v4);
        let ghost = [v0, v1, v2, v3].iter().position(|v| *v == Self::GHOST);

        let d0 = ghost.filter(|g| *g != 0)
            .map_or_else(|| (p2 - p3).cross(&(p1 - p3)).normalize().dot(&(p4 - p3)), |_| -f64::INFINITY);
        let d1 = ghost.filter(|g| *g != 1)
            .map_or_else(|| (p3 - p2).cross(&(p0 - p2)).normalize().dot(&(p4 - p2)), |_| -f64::INFINITY);
        let d2 = ghost.filter(|g| *g != 2)
            .map_or_else(|| (p0 - p1).cross(&(p3 - p1)).normalize().dot(&(p4 - p1)), |_| -f64::INFINITY);
        let d3 = ghost.filter(|g| *g != 3)
            .map_or_else(|| (p1 - p0).cross(&(p2 - p0)).normalize().dot(&(p4 - p0)), |_| -f64::INFINITY);
        d0.max(d1).max(d2).max(d3) < self.tolerance
    }

    /// Gets whether the triangle formed by the first 3 points intersects the edge formed by the last 2 points.
    /// Any vertex can be the ghost vertex.
    pub fn tri_intersects_edge_sim(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
        let positive = self.orient_3d_sim(v2, v1, v0, v3);
        self.orient_3d_sim(v0, v1, v2, v4) == positive &&
        self.orient_3d_sim(v0, v1, v3, v4) == positive &&
        self.orient_3d_sim(v1, v2, v3, v4) == positive &&
        self.orient_3d_sim(v2, v0, v3, v4) == positive
    }

    /// Gets unnormalized barycentric coordinates of the intersection point of a triangle and an edge.
    /// Any vertex can be the ghost vertex.
    pub fn tri_edge_bary(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> [f64; 5] {
        let bary0 = self.orient_3d(v1, v2, v3, v4);
        let bary1 = self.orient_3d(v2, v0, v3, v4);
        let bary2 = self.orient_3d(v0, v1, v3, v4);
        let bary3 = self.orient_3d(v0, v1, v2, v4);
        let bary4 = self.orient_3d(v2, v1, v0, v3);
        [bary0, bary1, bary2, bary3, bary4]
    }

    /// Gets whether the triangle formed by the first 3 points strictly intersects the edge formed by the last 2 points.
    /// Any vertex can be the ghost vertex.
    pub fn tri_intersects_edge(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
        let bary = self.tri_edge_bary(v0, v1, v2, v3, v4);
        bary.iter().all(|c| *c > 0.0) || bary.iter().all(|c| *c < 0.0)
    }

    /// Gets unnormalized barycentric coordinates of the intersection point of an edge and an edge.
    /// Any vertex can be the ghost vertex.
    pub fn edge_edge_bary(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
    ) -> [f64; 4] {
        let p0 = self.idx_ori(v0);
        let p1 = self.idx_ori(v1);
        let p2 = self.idx_ori(v2);
        let p3 = self.idx_ori(v3);

        // Cross product to determine if the lines are parallel
        let cross = Vec3::new(
            robust_geo::cross_2d(p0.yz(), p1.yz(), p2.yz(), p3.yz()),
            robust_geo::cross_2d(p0.zx(), p1.zx(), p2.zx(), p3.zx()),
            robust_geo::cross_2d(p0.xy(), p1.xy(), p2.xy(), p3.xy()),
        );

        if cross == Vec3::from_element(0.0) {
            [0.0; 4]
        } else {
            // Reduce to edge-triangle intersection
            let p4 = (if v2 == Self::GHOST {p3} else {p2}) + cross;
            let bary0 = (if [v1, v2, v3].contains(&Self::GHOST) {-1.0} else {1.0}) * robust_geo::orient_3d(p2, p1, p3, p4);
            let bary1 = (if [v2, v0, v3].contains(&Self::GHOST) {-1.0} else {1.0}) * robust_geo::orient_3d(p3, p0, p2, p4);
            let bary2 = (if [v0, v1, v3].contains(&Self::GHOST) {-1.0} else {1.0}) * robust_geo::orient_3d(p1, p3, p0, p4);
            let bary3 = (if [v0, v1, v2].contains(&Self::GHOST) {-1.0} else {1.0}) * robust_geo::orient_3d(p0, p2, p1, p4);
            [bary0, bary1, bary2, bary3]
        }
    }

    /// Gets whether the edge formed by the first 2 points strictly intersects the edge formed by the last 2 points.
    /// Any vertex can be the ghost vertex.
    pub fn edge_intersects_edge(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
    ) -> bool {
        if self.orient_3d(v0, v1, v2, v3) != 0.0 {
            return false;
        }
        let bary = self.edge_edge_bary(v0, v1, v2, v3);
        bary.iter().all(|c| *c > 0.0) || bary.iter().all(|c| *c < 0.0)
    }

    /// Gets the intersections along some edge in order from v0 to v1.
    /// v0 and v1 both have to be solid vertices.
    pub fn edge_intersections_raw(&self, v0: VertexId, v1: VertexId) -> Vec<EdgeIntersectionRaw> {
        fn tri_raw<V, T, I: IntoIterator<Item = TetWalker>>(mesh: &TetMesh<V, T>, walkers: I, v0: VertexId, v1: VertexId)
        -> (Option<(TetWalker, EdgeIntersectionRaw)>, bool) {
            let mut end = false;

            (walkers.into_iter().map(|walker| {
                let tri = walker.tri(mesh);
                if tri.contains(&v1) || end {
                    end = true;
                    return (walker, [0.0; 5])
                }

                let result = (walker, mesh.tri_edge_bary(tri[0], tri[1], tri[2], v0, v1));
                //println!("tri {}-{}-{}, v0 {}, v1 {}, bary {:?}", tri[0], tri[1], tri[2], v0, v1, result.1);
                result
            }).collect::<Vec<_>>().into_iter().find(|(_, bary)| // Second pass
                !end && bary.iter().all(|coord| *coord > 0.0) || bary.iter().all(|coord| *coord < 0.0))
            .map(|(walker, bary)| (
                walker,
                EdgeIntersectionRaw { other: EdgeOrTri::Tri(walker.tri(mesh)), bary: [bary[3], bary[4]]}
            )), end)
        }

        fn edge_raw<V, T, I: IntoIterator<Item = TetWalker>>(mesh: &TetMesh<V, T>, walkers: I, v0: VertexId, v1: VertexId)
        -> Option<(TetWalker, EdgeIntersectionRaw)> {
            walkers.into_iter().filter(|walker| {
                // Coplanar edges only
                let edge = walker.edge(mesh);
                let result = mesh.orient_3d(edge[0], edge[1], v0, v1);
                //println!("edge {}-{}, v0 {}, v1 {}, coplanarity test {}", edge[0], edge[1], v0, v1, result);
                result == 0.0
            }).map(|walker| {
                // No strict triangle intersections; look for exact edge intersection
                let edge = walker.edge(mesh);
                let result = (walker, mesh.edge_edge_bary(edge[0], edge[1], v0, v1));
                //println!("edge {}-{}, v0 {}, v1 {}, bary {:?}", edge[0], edge[1], v0, v1, result.1);
                result
            }).find(|(_, bary)|
                bary.iter().all(|coord| *coord > 0.0) || bary.iter().all(|coord| *coord < 0.0))
            .map(|(walker, bary)| (
                walker,
                EdgeIntersectionRaw { other: EdgeOrTri::Edge(walker.edge(mesh)), bary: [bary[2], bary[3]]}
            ))
        }

        // To avoid going backwards when searching from an edge.
        let mut prev = EdgeOrTri::Edge([VertexId::invalid(); 2]);

        // Look among triangles for strict intersection, then among edges
        let (mut first, end) = tri_raw(self, 
            self.walkers_from_vertex(v0).map(|w| w.to_perm(Permutation::_3210)), v0, v1);
        if !end {
            first = first.or_else(|| edge_raw(self, self.walkers_from_vertex(v0).flat_map(|walker| vec![
                walker.to_perm(Permutation::_3210),
                walker.to_perm(Permutation::_2130),
                walker.to_perm(Permutation::_1320),
            ]), v0, v1));
        }

        iter::successors(first, |(walker, inter)| {
            let result = match inter.other {
                EdgeOrTri::Tri(_) => {
                    let walker = walker.to_adj_ae(self);

                    // Check all 3 other triangles, then all 3 other edges.
                    let (res, end) = tri_raw(self, vec![walker.to_twin_edge(), walker.to_opp_edge(), walker.to_perm(Permutation::_3210)],
                        v0, v1);
                    if end {
                        None
                    } else {
                        res.or_else(|| edge_raw(self, vec![
                            walker.to_perm(Permutation::_0312),
                            walker.to_perm(Permutation::_1320),
                            walker.to_perm(Permutation::_2301),
                        ], v0, v1))
                    }
                },

                EdgeOrTri::Edge(_) => {
                    // Check all triangles on the boundary of the edge shell, then all other edges of the shell.

                    let (res, end) = tri_raw(self, walker.edge_ring(self)
                        .flat_map(|w| vec![w.to_opp_edge(), w.to_perm(Permutation::_3210)])
                        .filter(|w| Some(sorted_3(w.tri(self))) != prev.as_tri()), v0, v1);
                    
                    if end {
                        None
                    } else {
                        res.or_else(|| edge_raw(self, walker.edge_ring(self).flat_map(|w| vec![
                                w.to_perm(Permutation::_0312),
                                w.to_perm(Permutation::_1320),
                                w.to_perm(Permutation::_2301),
                            ].into_iter().filter(|w| Some(sorted_2(w.edge(self))) != prev.as_edge())), v0, v1))
                    }
                },
            };
            prev = match inter.other {
                EdgeOrTri::Edge(edge) => EdgeOrTri::Edge(sorted_2(edge)),
                EdgeOrTri::Tri(tri) => EdgeOrTri::Tri(sorted_3(tri)),
            };
            result
        }).map(|(_, inter)| inter).collect()
    }
    
    /// Gets the intersections along some edge in order from v0 to v1.
    /// v0 and v1 both have to be solid vertices.
    pub fn edge_intersections(&self, v0: VertexId, v1: VertexId) -> Vec<EdgeIntersection> {
        let p0 = self[v0].position();
        let p1 = self[v1].position();

        self.edge_intersections_raw(v0, v1).into_iter().map(|inter| {
            let sum = inter.bary[0] + inter.bary[1];
            EdgeIntersection {
                other: inter.other,
                point: p0 + (p1 - p0) * inter.bary[1] / sum
            }
        }).collect::<Vec<_>>()
    }

    /// Gets the intersection along some edge closest to the edge's midpoint.
    /// v0 and v1 both have to be solid vertices.
    pub fn edge_mid_intersection(&self, v0: VertexId, v1: VertexId) -> Option<EdgeIntersection> {
        let p0 = self[v0].position();
        let p1 = self[v1].position();

        self.edge_intersections_raw(v0, v1).into_iter().min_by_key(|inter| {
            FloatOrd((inter.bary[0] - inter.bary[1]).abs())
        }).map(|inter| {
            let sum = inter.bary[0] + inter.bary[1];
            EdgeIntersection {
                other: inter.other,
                point: p0 + (p1 - p0) * inter.bary[1] / sum
            }
        })
    }

    /// Gets whether the last point is in the tet first 4 points.
    /// The tet must be oriented positive.
    pub fn tet_intersects_vertex(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
        self.orient_3d_sim(v0, v1, v2, v4) &&
        self.orient_3d_sim(v3, v2, v1, v4) &&
        self.orient_3d_sim(v2, v3, v0, v4) &&
        self.orient_3d_sim(v1, v0, v3, v4)
    }

    /// Gets the intersections along multiple triangles.
    /// The triangle vertices have to be solid.
    pub fn tris_intersections_raw(&self, tris: &[[VertexId; 3]], start_edge: [VertexId; 2]) -> Vec<(TriIntersectionRaw, [VertexId; 3])> {
        fn edge_raw<V, T, I: IntoIterator<Item = TetWalker>>(mesh: &TetMesh<V, T>, walkers: I, tris: &[[VertexId; 3]])
        -> Option<(TetWalker, TriIntersectionRaw, [VertexId; 3])> {
            walkers.into_iter().flat_map(|walker| {
                // Find first triangle that intersects
                tris.iter().map(move |[v0, v1, v2]| {
                    let edge = walker.edge(mesh);
                    let result = (walker, mesh.tri_edge_bary(*v0, *v1, *v2, edge[0], edge[1]), [*v0, *v1, *v2]);
                    //println!("edge {}-{}, v0 {}, v1 {}, v2 {}, bary {:?}", edge[0], edge[1], v0, v1, v2, result.1);
                    result
                })
            }).find(|(walker, bary, _)| {
                // Allow zeros to account for triangulation edges, but make sure vertices of tris don't count as intersections

                let edge = walker.edge(mesh);
                !tris.iter().any(|tri| tri.contains(&edge[0]) || tri.contains(&edge[1])) &&
                    bary.iter().all(|coord| *coord >= 0.0) || bary.iter().all(|coord| *coord <= 0.0)
            }).map(|(walker, bary, tri)| (
                walker,
                TriIntersectionRaw { edge: walker.edge(mesh), bary: [bary[0], bary[1], bary[2]]},
                tri,
            ))
        }

        let walker = self.walker_from_edge(start_edge).unwrap();
        
        // Find first intersection
        if let Some((walker, inter, tri)) = edge_raw(self, walker.edge_ring(self).map(TetWalker::to_opp_edge),
            tris
        ) {
            let mut visited = FnvHashSet::default();
            visited.insert(sorted_2(walker.edge(self)));

            let mut to_search = walker.edge_ring(self).flat_map(|w|
                vec![w.to_nfe(), w.to_pfe()]).collect::<Vec<_>>();

            let mut inters = vec![(inter, tri)];

            while let Some(walker) = to_search.pop() {
                if visited.insert(sorted_2(walker.edge(self))) {
                    if let Some((walker, inter, tri)) = edge_raw(self, iter::once(walker), tris) {
                        inters.push((inter, tri));
                        to_search.extend(walker.edge_ring(self).flat_map(|w|
                            vec![w.to_nfe(), w.to_pfe()]));
                    }
                }
            }

            inters
        } else {
            vec![]
        }
    }

    /// Gets the intersections along some triangle.
    /// v0, v1, v2 have to be solid vertices.
    /// Assumes that all edges of the triangle are present.
    pub fn tri_intersections_raw(&self, v0: VertexId, v1: VertexId, v2: VertexId) -> Vec<TriIntersectionRaw> {
        fn edge_raw<V, T, I: IntoIterator<Item = TetWalker>>(mesh: &TetMesh<V, T>, walkers: I, v0: VertexId, v1: VertexId, v2: VertexId)
        -> Option<(TetWalker, TriIntersectionRaw)> {
            walkers.into_iter().map(|walker| {
                let edge = walker.edge(mesh);
                let result = (walker, mesh.tri_edge_bary(v0, v1, v2, edge[0], edge[1]));
                //println!("edge {}-{}, v0 {}, v1 {}, v2 {}, bary {:?}", edge[0], edge[1], v0, v1, v2, result.1);
                result
            }).find(|(_, bary)|
                bary.iter().all(|coord| *coord > 0.0) || bary.iter().all(|coord| *coord < 0.0))
            .map(|(walker, bary)| (
                walker,
                TriIntersectionRaw { edge: walker.edge(mesh), bary: [bary[0], bary[1], bary[2]]}
            ))
        }

        let walker = self.walker_from_edge([v0, v1]).unwrap();
        
        // Find first intersection
        if let Some((walker, inter)) = edge_raw(self, walker.edge_ring(self).map(TetWalker::to_opp_edge),
            v0, v1, v2
        ) {
            let mut visited = FnvHashSet::default();
            visited.insert(sorted_2(walker.edge(self)));

            let mut to_search = walker.edge_ring(self).flat_map(|w|
                vec![w.to_nfe(), w.to_pfe()]).collect::<Vec<_>>();

            let mut inters = vec![inter];

            while let Some(walker) = to_search.pop() {
                if visited.insert(sorted_2(walker.edge(self))) {
                    if let Some((walker, inter)) = edge_raw(self, iter::once(walker), v0, v1, v2) {
                        inters.push(inter);
                        to_search.extend(walker.edge_ring(self).flat_map(|w|
                            vec![w.to_nfe(), w.to_pfe()]));
                    }
                }
            }

            inters
        } else {
            vec![]
        }
    }
    
    /// Gets the intersections along some triangles.
    pub fn tris_intersections(&self, tris: &[[VertexId; 3]], start_edge: [VertexId; 2]) -> Vec<TriIntersection> {
        self.tris_intersections_raw(tris, start_edge).into_iter().map(|(inter, [v0, v1, v2])| {
            let p0 = self[v0].position();
            let p1 = self[v1].position();
            let p2 = self[v2].position();

            let sum = inter.bary[0] + inter.bary[1] + inter.bary[2];
            TriIntersection {
                edge: inter.edge,
                point: p0 + (p1 - p0) * inter.bary[1] / sum + (p2 - p0) * inter.bary[2] / sum
            }
        }).collect::<Vec<_>>()
    }
    
    /// Gets the intersections along some triangle.
    /// v0, v1, v2 have to be solid vertices.
    /// Assumes that all edges of the triangle are present.
    pub fn tri_intersections(&self, v0: VertexId, v1: VertexId, v2: VertexId) -> Vec<TriIntersection> {
        let p0 = self[v0].position();
        let p1 = self[v1].position();
        let p2 = self[v2].position();

        self.tri_intersections_raw(v0, v1, v2).into_iter().map(|inter| {
            let sum = inter.bary[0] + inter.bary[1] + inter.bary[2];
            TriIntersection {
                edge: inter.edge,
                point: p0 + (p1 - p0) * inter.bary[1] / sum + (p2 - p0) * inter.bary[2] / sum
            }
        }).collect::<Vec<_>>()
    }

    /// Gets the amount of time that `center` has to travel toward `target` to not be star-shaped relative to
    /// `ball_vertex`.
    /// Returns None if no such time exists.
    pub fn time_until_not_star_shaped(&self, ball_vertex: VertexId, center: VertexId, target: VertexId) -> Option<f64> {
        self.walkers_from_vertex(ball_vertex).map(|walker| {
            let tri = walker.opp_tri(self);
            let bary = self.tri_edge_bary(tri[0], tri[1], tri[2], center, target);
            bary
        }).filter(|bary| bary.iter().all(|c| *c >= 0.0) || bary.iter().all(|c| *c <= 0.0))
            .map(|bary| bary[4] / (bary[3] + bary[4]))
            .filter(|time| *time > 0.0 && time.is_finite())
            .min_by_key(|time| FloatOrd(*time))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{AbsDiffEq, UlpsEq, assert_ulps_eq};

    #[derive(Clone, Debug, PartialEq)]
    struct EI(Vec<EdgeIntersection>);

    impl AbsDiffEq<EI> for EI {
        type Epsilon = <f64 as AbsDiffEq<f64>>::Epsilon;

        fn default_epsilon() -> Self::Epsilon {
            f64::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &EI, epsilon: Self::Epsilon) -> bool {
            self.0.iter().zip(other.0.iter()).all(|(a, b)| a.point.abs_diff_eq(&b.point, epsilon))
        }
    }

    impl UlpsEq<EI> for EI {
        fn default_max_ulps() -> u32 {
            f64::default_max_ulps()
        }

        fn ulps_eq(&self, other: &EI, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
            self.0.iter().zip(other.0.iter()).all(|(a, b)| a.point.ulps_eq(&b.point, epsilon, max_ulps))
        }
    }
    
    fn v(n: IdType) -> VertexId {
        VertexId(n)
    }

    #[derive(Clone, Debug, PartialEq)]
    struct TI(Vec<TriIntersection>);

    impl AbsDiffEq<TI> for TI {
        type Epsilon = <f64 as AbsDiffEq<f64>>::Epsilon;

        fn default_epsilon() -> Self::Epsilon {
            f64::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &TI, epsilon: Self::Epsilon) -> bool {
            self.0.iter().zip(other.0.iter()).all(|(a, b)| a.point.abs_diff_eq(&b.point, epsilon))
        }
    }

    impl UlpsEq<TI> for TI {
        fn default_max_ulps() -> u32 {
            f64::default_max_ulps()
        }

        fn ulps_eq(&self, other: &TI, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
            self.0.iter().zip(other.0.iter()).all(|(a, b)| a.point.ulps_eq(&b.point, epsilon, max_ulps))
        }
    }

    #[test]
    fn test_edge_intersections_none() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 0.2), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(3), v(4));
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![]));
    }

    #[test]
    fn test_edge_intersections_one_tri() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 2.0), ()),
                (Pt3::new(0.0, 0.0, -2.0), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(3), v(4));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Tri([v(0), v(1), v(2)]), Pt3::new(0.0, 0.0, 0.0)),
        ]));
    }

    #[test]
    fn test_edge_intersections_one_edge() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(3.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 2.0, 0.0), ()),
                (Pt3::new(-3.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, -2.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(0.0, 0.0, -1.0), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(3), v(1));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Edge([v(4), v(5)]), Pt3::new(0.0, 0.0, 0.0)),
        ]));
    }

    #[test]
    fn test_edge_intersections_multi_edge() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 3.0, 0.0), ()),
                (Pt3::new(2.0, 1.0, 0.0), ()),
                (Pt3::new(2.0, -1.0, 0.0), ()),
                (Pt3::new(1.0, -3.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(1), v(4));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Edge([v(0), v(2)]), Pt3::new(1.0, 0.5, 0.0)),
            EdgeIntersection::new(EdgeOrTri::Edge([v(0), v(3)]), Pt3::new(1.0, -0.5, 0.0)),
        ]));
    }

    #[test]
    fn test_edge_intersections_edge_tri() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 3.0, 0.0), ()),
                (Pt3::new(2.0, 1.0, 0.0), ()),
                (Pt3::new(2.0, -1.0, 0.1), ()),
                (Pt3::new(1.0, -3.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(1), v(4));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Edge([v(0), v(2)]), Pt3::new(1.0, 0.5, 0.0)),
            EdgeIntersection::new(EdgeOrTri::Tri([v(0), v(3), v(5)]), Pt3::new(1.0, -0.5, 0.0)),
        ]), max_ulps = 8);
    }

    #[test]
    fn test_edge_intersections_tri_edge() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 3.0, 0.0), ()),
                (Pt3::new(2.0, 1.0, 0.1), ()),
                (Pt3::new(2.0, -1.0, 0.0), ()),
                (Pt3::new(1.0, -3.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(1), v(4));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Tri([v(0), v(2), v(5)]), Pt3::new(1.0, 0.5, 0.0)),
            EdgeIntersection::new(EdgeOrTri::Edge([v(0), v(3)]), Pt3::new(1.0, -0.5, 0.0)),
        ]));
    }

    #[test]
    fn test_edge_intersections_multi_tri() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 3.0, 0.0), ()),
                (Pt3::new(2.0, 1.0, 0.1), ()),
                (Pt3::new(2.0, -1.0, 0.1), ()),
                (Pt3::new(1.0, -3.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let inters = mesh.edge_intersections(v(1), v(4));
        #[rustfmt::skip]
        assert_ulps_eq!(EI(inters.into_iter().map(EdgeIntersection::sorted).collect::<Vec<_>>()), EI(vec![
            EdgeIntersection::new(EdgeOrTri::Tri([v(0), v(2), v(5)]), Pt3::new(1.0, 0.5, 0.0)),
            EdgeIntersection::new(EdgeOrTri::Tri([v(0), v(3), v(5)]), Pt3::new(1.0, -0.5, 0.0)),
        ]), max_ulps = 8);
    }

    #[test]
    fn test_edges_intersecting_tri_none() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 2.0), ()),
                (Pt3::new(0.0, 0.0, -2.0), ()),
            ],
            0.0,
            || (),
        );

        let edges = mesh.tri_intersections(v(0), v(1), v(2));
        assert_ulps_eq!(TI(edges.into_iter().map(TriIntersection::sorted).collect::<Vec<_>>()), TI(vec![]));
    }

    #[test]
    fn test_edges_intersecting_tri_one() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 0.2), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let edges = mesh.tri_intersections(v(0), v(1), v(2));
        assert_ulps_eq!(TI(edges.into_iter().map(TriIntersection::sorted).collect::<Vec<_>>()), TI(vec![
            TriIntersection::new([v(3), v(4)], Pt3::new(0.0, 0.0, 0.0))
        ]));
    }

    #[test]
    fn test_edges_intersecting_tri_some() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(-0.1, 0.0, 0.2), ()),
                (Pt3::new(0.1, 0.0, 0.2), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );

        let edges = mesh.tri_intersections(v(0), v(1), v(2));
        assert_ulps_eq!(TI(edges.into_iter().map(TriIntersection::sorted).collect::<Vec<_>>()), TI(vec![
            // Order shouldn't matter, but I want to move on to Steiner point suppression already.
            TriIntersection::new([v(4), v(5)], Pt3::new(0.05, 0.0, 0.0)),
            TriIntersection::new([v(3), v(5)], Pt3::new(-0.05, 0.0, 0.0)),
        ]));
    }
}