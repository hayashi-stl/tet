use std::hash::{Hash, Hasher};
use std::ops::Index;

use crate::id_map::{IdMap, IdType};
use crate::{Pt3, TetId, VertexId};
use fnv::FnvHashSet;
use simplicity as sim;

/// A vertex stores a tet that it's part of and its position.
#[derive(Clone, Debug)]
pub struct Vertex<V> {
    tet: TetWalker,
    position: Pt3,
    value: V,
}

impl<V> Vertex<V> {
    fn new(tet: TetWalker, position: Pt3, value: V) -> Self {
        Self {
            tet,
            position,
            value,
        }
    }

    /// Gets the position of this vertex.
    pub fn position(&self) -> Pt3 {
        self.position
    }

    /// Gets the value of this vertex.
    pub fn value(&self) -> &V {
        &self.value
    }
}

/// A tet stores the vertices that it's part of (in the correct orientation)
/// and the opposite tets of each vertex.
///
/// At most 1 vertex id is allowed to be invalid. If there is an
/// invalid vertex id, then this is a ghost tet.
///
///
/// Walkers store the index of the edge of the opposite triangle
/// that is the twin of the first edge of the current triangle.
/// This idea is borrowed from TetGen.
#[derive(Clone, Debug)]
pub struct Tet<T> {
    vertices: [VertexId; 4],
    opp_tets: [TetWalker; 4],
    value: T,
}

impl<T> Tet<T> {
    fn new(vertices: [VertexId; 4], opp_tets: [TetWalker; 4], value: T) -> Self {
        Self {
            vertices,
            opp_tets,
            value,
        }
    }

    /// Gets the vertices of this tet. The ghost vertex may be included.
    pub fn vertices(&self) -> [VertexId; 4] {
        self.vertices
    }

    /// Gets the value of this tet.
    pub fn value(&self) -> &T {
        &self.value
    }
}

/// A boundary triangle processed by `TetWalker::boundary_and_enclosed`.
/// Use `self.0` to access the tet walker.
#[derive(Clone, Copy, Debug)]
struct BoundaryTri(TetWalker);

impl PartialEq for BoundaryTri {
    fn eq(&self, other: &BoundaryTri) -> bool {
        self.0.id() == other.0.id() && self.0.edge % 4 == other.0.edge % 4
    }
}

impl Eq for BoundaryTri {}

impl Hash for BoundaryTri {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.id().hash(state);
        (self.0.edge % 4).hash(state);
    }
}

const fn permutation_arr(inv: bool) -> [[u32; 12]; 12] {
    let mut arr = [[0; 12]; 12];
    let perms = [
        [3, 2, 1, 0],
        [2, 3, 0, 1],
        [1, 0, 3, 2],
        [0, 1, 2, 3],
        [2, 1, 3, 0],
        [3, 0, 2, 1],
        [0, 3, 1, 2],
        [1, 2, 0, 3],
        [1, 3, 2, 0],
        [0, 2, 3, 1],
        [3, 1, 0, 2],
        [2, 0, 1, 3],
    ];

    let mut p = 0;
    while p < 12 {
        let mut i = 0;
        while i < 12 {
            let mut perm = perms[p];
            if inv {
                let mut inv = [0; 4];
                inv[perm[0]] = 0;
                inv[perm[1]] = 1;
                inv[perm[2]] = 2;
                inv[perm[3]] = 3;
                perm = inv;
            }

            let edge = [perm[perms[i][0]], perm[perms[i][1]]];
            let mut new = 0;
            while new < 12 {
                if edge[0] == perms[new][0] && edge[1] == perms[new][1] {
                    arr[p][i] = new as u32;
                }
                new += 1;
            }

            i += 1;
        }
        p += 1;
    }

    arr
}

/// A permutation of the vertices of a tet walker.
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub enum Permutation {
    _3210,
    _2301,
    _1032,
    _0123,
    _2130,
    _3021,
    _0312,
    _1203,
    _1320,
    _0231,
    _3102,
    _2013,
}

/// A walker over the tetrahedrons of a tet mesh.
/// The edge index is a number in 0..12 indexing into
/// the edges of a tet.
///
/// The enumeration of the edges is borrowed from TetGen.
/// One edge of each triangle is enumerated,
/// then the next edge of each triangle respectively,
/// then the last edge of each triangle respectively.
///
/// So if the vertices of the tet are 0, 1, 2, 3, then the edges (and opposite edges) are
/// ```notrust
///  0: 3-2 (1-0)
///  1: 2-3 (0-1)
///  2: 1-0 (3-2)
///  3: 0-1 (2-3)
///  4: 2-1 (3-0)
///  5: 3-0 (2-1)
///  6: 0-3 (1-2)
///  7: 1-2 (0-3)
///  8: 1-3 (2-0)
///  9: 0-2 (3-1)
/// 10: 3-1 (0-2)
/// 11: 2-0 (1-3)
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct TetWalker {
    tet: TetId,
    edge: u32,
}

impl TetWalker {
    /// Index by permutation first.
    /// Permutation abcd takes 0123 → abcd
    const PERMUTATION_FWD: [[u32; 12]; 12] = permutation_arr(false);

    /// Index by permutation first.
    /// Permutation-inverse abcd takes abcd → 0123
    const PERMUTATION_INV: [[u32; 12]; 12] = permutation_arr(true);

    fn new(tet: TetId, edge_index: u32) -> Self {
        Self {
            tet,
            edge: edge_index,
        }
    }

    /// Gets the first vertex of the tet walker. This is the current vertex.
    pub fn first<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2][self.edge as usize]]
    }

    /// Gets the second vertex of the tet walker. This is the current edge's target.
    pub fn second<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[2, 3, 0, 1, 1, 0, 3, 2, 3, 2, 1, 0][self.edge as usize]]
    }

    /// Gets the third vertex of the tet walker.
    pub fn third<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1][self.edge as usize]]
    }

    /// Gets the fourth vertex of the tet walker. This is the vertex opposite the current triangle.
    pub fn fourth<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3][self.edge as usize]]
    }

    /// Gets the current edge of the tet walker.
    pub fn edge<V, T>(self, mesh: &Tets<V, T>) -> [VertexId; 2] {
        [self.first(mesh), self.second(mesh)]
    }

    /// Gets the opposite edge of the current edge on the tet walker's tet.
    pub fn opp_edge<V, T>(self, mesh: &Tets<V, T>) -> [VertexId; 2] {
        [self.third(mesh), self.fourth(mesh)]
    }

    /// Gets the current triangle of the tet walker.
    pub fn tri<V, T>(self, mesh: &Tets<V, T>) -> [VertexId; 3] {
        [self.first(mesh), self.second(mesh), self.third(mesh)]
    }

    /// Gets the opposite triangle of the current vertex on the tet walker's tet.
    /// The vertices are returned in [fourth, third, second] order.
    pub fn opp_tri<V, T>(self, mesh: &Tets<V, T>) -> [VertexId; 3] {
        [self.fourth(mesh), self.third(mesh), self.second(mesh)]
    }

    /// Gets the vertices of the current tet, respecting the orientation.
    pub fn tet<V, T>(self, mesh: &Tets<V, T>) -> [VertexId; 4] {
        [
            self.first(mesh),
            self.second(mesh),
            self.third(mesh),
            self.fourth(mesh),
        ]
    }

    /// Gets the current tet id.
    pub fn id(self) -> TetId {
        self.tet
    }

    crate::alias! {
        to_nfe,
        /// Advance to the next edge in the current triangle.
        pub fn to_next_tri_edge((mut) self: Self) -> Self {
            self.edge = (self.edge + 4) % 12;
            self
        }
    }

    crate::alias! {
        to_pfe,
        /// Advance to the previous edge in the current triangle.
        pub fn to_prev_tri_edge((mut) self: Self) -> Self {
            self.edge = (self.edge + 8) % 12;
            self
        }
    }

    crate::alias! {
        to_nve,
        /// Advance to the next edge from the current vertex. This is a counterclockwise rotation.
        pub fn to_next_vertex_edge((mut) self: Self) -> Self {
            self.edge = [10, 11, 8, 9, 1, 0, 3, 2, 7, 6, 5, 4][self.edge as usize];
            self
        }
    }

    crate::alias! {
        to_pve,
        /// Advance to the previous edge from the current vertex. This is a clockwise rotation.
        pub fn to_prev_vertex_edge((mut) self: Self) -> Self {
            self.edge = [5, 4, 7, 6, 11, 10, 9, 8, 2, 3, 0, 1][self.edge as usize];
            self
        }
    }

    /// Flip the current edge, keeping the current tet the same.
    pub fn to_twin_edge(mut self) -> Self {
        self.edge = [1, 0, 3, 2, 7, 6, 5, 4, 10, 11, 8, 9][self.edge as usize];
        self
    }

    /// Sets the current edge to the opposite edge of the current tet.
    pub fn to_opp_edge(mut self) -> Self {
        self.edge = [2, 3, 0, 1, 5, 4, 7, 6, 11, 10, 9, 8][self.edge as usize];
        self
    }

    /// Permutes the current vertex order with some forward permutation.
    /// Use when none of the same-tet operations are what you want.
    ///
    /// Special cases:
    ///
    /// Perm | Function
    /// -----|---------
    /// 0123 | (identity)
    /// 1203 | `to_next_tri_edge` aka `to_nfe`
    /// 2013 | `to_prev_tri_edge` aka `to_pfe`
    /// 0231 | `to_next_vertex_edge` aka `to_nve`
    /// 0312 | `to_prev_vertex_edge` aka `to_pve`
    /// 1032 | `to_twin_edge`
    /// 2301 | `to_opp_edge`
    pub fn to_perm(mut self, perm: Permutation) -> Self {
        // Don't question the indexing order here, it just works™
        self.edge = Self::PERMUTATION_FWD[self.edge as usize][perm as usize];
        self
    }

    crate::alias! {
        to_adj,
        /// Flip the current triangle, moving to the adjacent tet on that triangle.
        /// This flips the current edge.
        pub fn to_twin_tri<V, T>((mut) self: Self, mesh: &Tets<V, T>) -> Self {
            let div = self.edge / 4;
            self = mesh[self.id()].opp_tets[self.edge as usize % 4];
            match div {
                0 => self,
                1 => self.to_prev_tri_edge(),
                2 => self.to_next_tri_edge(),
                _ => unreachable!(),
            }
        }
    }

    crate::alias! {
        to_adj_ae,
        /// Flip the current triangle, moving to the adjacent tet on that triangle.
        /// This moves to some edge on that triangle.
        /// Use this method if you don't care which edge the walker ends up on.
        pub fn to_twin_tri_any_edge<V, T>(self: Self, mesh: &Tets<V, T>) -> Self {
            mesh[self.id()].opp_tets[self.edge as usize % 4]
        }
    }

    /// Sets the orientation such that the vertices returned by `self.vertices()`
    /// are in the same order as those returned by `mesh[self.id()].vertices()`.
    pub fn to_canon_tet(mut self) -> Self {
        self.edge = 3;
        self
    }

    /// Performs a 1-to-4 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Returns walkers of the 4 tets added. Each walker's 4th vertex is the new vertex.
    pub fn flip14_unchecked<V: Clone, T: Clone>(
        self,
        mesh: &mut Tets<V, T>,
        vertex: VertexId,
    ) -> [TetWalker; 4] {
        let id0 = self.id();
        let id1 = mesh.tets.insert(mesh.tets[id0].clone());
        let id2 = mesh.tets.insert(mesh.tets[id0].clone());
        let id3 = mesh.tets.insert(mesh.tets[id0].clone());

        // The vertex didn't have a tet before.
        mesh.vertices[vertex].tet = self;

        // Update tets
        // Note that adjacent tets need to change their opposite tet ids,
        // but not their opposite tet edge indexes because they refer to the same face.
        mesh.tets[id0].vertices[0] = vertex;
        mesh.tets[id0].opp_tets[1] = TetWalker::new(id1, 0);
        mesh.tets[id0].opp_tets[2] = TetWalker::new(id2, 4);
        mesh.tets[id0].opp_tets[3] = TetWalker::new(id3, 8);

        mesh.tets[id1].vertices[1] = vertex;
        mesh.tets[id1].opp_tets[2] = TetWalker::new(id2, 9);
        mesh.tets[id1].opp_tets[3] = TetWalker::new(id3, 5);
        mesh.tets[id1].opp_tets[0] = TetWalker::new(id0, 1);

        mesh.tets[id2].vertices[2] = vertex;
        mesh.tets[id2].opp_tets[3] = TetWalker::new(id3, 2);
        mesh.tets[id2].opp_tets[0] = TetWalker::new(id0, 6);
        mesh.tets[id2].opp_tets[1] = TetWalker::new(id1, 10);

        mesh.tets[id3].vertices[3] = vertex;
        mesh.tets[id3].opp_tets[0] = TetWalker::new(id0, 11);
        mesh.tets[id3].opp_tets[1] = TetWalker::new(id1, 7);
        mesh.tets[id3].opp_tets[2] = TetWalker::new(id2, 3);

        // Update adjacent tets' opposite indexes
        let ids = [id0, id1, id2, id3];
        for i in 1..4 {
            let walker = mesh[ids[i]].opp_tets[i];
            mesh.tets[walker.id()].opp_tets[walker.edge as usize % 4].tet = ids[i];
        }

        [
            TetWalker::new(id0, 0),
            TetWalker::new(id1, 1),
            TetWalker::new(id2, 2),
            TetWalker::new(id3, 3),
        ]
    }

    /// Performs a 2-to-3 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Returns walkers of the 3 tets created. Each walker's opposite edge is the new edge.
    pub fn flip23_unchecked<V: Clone, T: Clone>(self, mesh: &mut Tets<V, T>) -> [TetWalker; 3] {
        // Walker on other tri
        let other = self.to_adj(&mesh);

        let id0 = self.id();
        let id1 = other.id();
        let id2 = mesh.tets.insert(mesh.tets[id0].clone());

        // Fix tet walkers on adjacent tets
        let walk0t = self.to_perm(Permutation::_1032);
        let adj0t = walk0t.to_adj(&mesh);
        let walker = &mut mesh.tets[adj0t.id()].opp_tets[adj0t.edge as usize % 4];
        *walker = TetWalker::new(
            id0,
            Self::PERMUTATION_INV[walk0t.edge as usize][walker.edge as usize],
        );

        let walk1t = self.to_perm(Permutation::_0231);
        let adj1t = walk1t.to_adj(&mesh);
        let walker = &mut mesh.tets[adj1t.id()].opp_tets[adj1t.edge as usize % 4];
        *walker = TetWalker::new(
            id1,
            Self::PERMUTATION_INV[walk1t.edge as usize][walker.edge as usize],
        );

        let walk2t = self.to_perm(Permutation::_2130);
        let adj2t = walk2t.to_adj(&mesh);
        let walker = &mut mesh.tets[adj2t.id()].opp_tets[adj2t.edge as usize % 4];
        *walker = TetWalker::new(
            id2,
            Self::PERMUTATION_INV[walk2t.edge as usize][walker.edge as usize],
        );

        // Bottom tets; permutations are a little different
        let adj0b = other.to_perm(Permutation::_1032).to_adj(&mesh);
        let walker = &mut mesh.tets[adj0b.id()].opp_tets[adj0b.edge as usize % 4];
        *walker = TetWalker::new(
            id0,
            Self::PERMUTATION_INV[other.edge as usize][walker.edge as usize],
        );

        let adj1b = other.to_perm(Permutation::_2130).to_adj(&mesh);
        let walker = &mut mesh.tets[adj1b.id()].opp_tets[adj1b.edge as usize % 4];
        *walker = TetWalker::new(
            id1,
            Self::PERMUTATION_INV[other.to_nfe().edge as usize][walker.edge as usize],
        );

        let adj2b = other.to_perm(Permutation::_0231).to_adj(&mesh);
        let walker = &mut mesh.tets[adj2b.id()].opp_tets[adj2b.edge as usize % 4];
        *walker = TetWalker::new(
            id2,
            Self::PERMUTATION_INV[other.to_pfe().edge as usize][walker.edge as usize],
        );

        // Introduce the 3 new tets
        let ([v1, v0, v2, v3], v4) = (self.tet(&mesh), other.fourth(&mesh));

        let tet = &mut mesh.tets[id0];
        tet.vertices = [v0, v1, v3, v4];
        tet.opp_tets[0] = TetWalker::new(id1, 1);
        tet.opp_tets[1] = TetWalker::new(id2, 0);
        tet.opp_tets[2] = TetWalker::new(adj0b.id(), adj0b.edge);
        tet.opp_tets[3] = TetWalker::new(adj0t.id(), adj0t.edge);

        let tet = &mut mesh.tets[id1];
        tet.vertices = [v1, v2, v3, v4];
        tet.opp_tets[0] = TetWalker::new(id2, 1);
        tet.opp_tets[1] = TetWalker::new(id0, 0);
        tet.opp_tets[2] = TetWalker::new(adj1b.id(), adj1b.edge);
        tet.opp_tets[3] = TetWalker::new(adj1t.id(), adj1t.edge);

        let tet = &mut mesh.tets[id2];
        tet.vertices = [v2, v0, v3, v4];
        tet.opp_tets[0] = TetWalker::new(id0, 1);
        tet.opp_tets[1] = TetWalker::new(id1, 0);
        tet.opp_tets[2] = TetWalker::new(adj2b.id(), adj2b.edge);
        tet.opp_tets[3] = TetWalker::new(adj2t.id(), adj2t.edge);

        [
            TetWalker::new(id0, 3),
            TetWalker::new(id1, 3),
            TetWalker::new(id2, 3),
        ]
    }

    /// Performs a 3-to-2 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Returns walkers for the 2 tets created by the flip. Each walker's opposite triangle is the new triangle.
    pub fn flip32_unchecked<V: Clone, T: Clone>(self, mesh: &mut Tets<V, T>) -> [TetWalker; 2] {
        let other1 = self.to_twin_edge().to_adj(&mesh);
        let other2 = other1.to_twin_edge().to_adj(&mesh);

        let id0 = self.id();
        let id1 = other1.id();

        // Fix tet walkers on adjacent tets
        let adj0t = self.to_opp_edge().to_adj(&mesh);
        let walker = &mut mesh.tets[adj0t.id()].opp_tets[adj0t.edge as usize % 4];
        *walker = TetWalker::new(
            id0,
            Self::PERMUTATION_INV[self.to_perm(Permutation::_3210).edge as usize]
                [walker.edge as usize],
        );

        let adj1t = other1.to_opp_edge().to_adj(&mesh);
        let walker = &mut mesh.tets[adj1t.id()].opp_tets[adj1t.edge as usize % 4];
        *walker = TetWalker::new(
            id0,
            Self::PERMUTATION_INV[other1.to_perm(Permutation::_2130).edge as usize]
                [walker.edge as usize],
        );

        let adj2t = other2.to_opp_edge().to_adj(&mesh);
        let walker = &mut mesh.tets[adj2t.id()].opp_tets[adj2t.edge as usize % 4];
        *walker = TetWalker::new(
            id0,
            Self::PERMUTATION_INV[other2.to_perm(Permutation::_1320).edge as usize]
                [walker.edge as usize],
        );

        // Bottom
        let adj0b = self.to_perm(Permutation::_3210).to_adj(&mesh);
        let walker = &mut mesh.tets[adj0b.id()].opp_tets[adj0b.edge as usize % 4];
        *walker = TetWalker::new(
            id1,
            Self::PERMUTATION_INV[self.to_perm(Permutation::_2301).edge as usize]
                [walker.edge as usize],
        );

        let adj1b = other1.to_perm(Permutation::_3210).to_adj(&mesh);
        let walker = &mut mesh.tets[adj1b.id()].opp_tets[adj1b.edge as usize % 4];
        *walker = TetWalker::new(
            id1,
            Self::PERMUTATION_INV[other1.to_perm(Permutation::_0231).edge as usize]
                [walker.edge as usize],
        );

        let adj2b = other2.to_perm(Permutation::_3210).to_adj(&mesh);
        let walker = &mut mesh.tets[adj2b.id()].opp_tets[adj2b.edge as usize % 4];
        *walker = TetWalker::new(
            id1,
            Self::PERMUTATION_INV[other2.to_perm(Permutation::_3021).edge as usize]
                [walker.edge as usize],
        );

        // Introduce the 2 new tets
        let ([v3, v4, v1, v0], v2) = (self.tet(&mesh), other1.fourth(&mesh));

        let tet = &mut mesh.tets[id0];
        tet.vertices = [v0, v1, v2, v3];
        tet.opp_tets[0] = TetWalker::new(adj2t.id(), adj2t.to_nfe().edge);
        tet.opp_tets[1] = TetWalker::new(adj1t.id(), adj1t.to_pfe().edge);
        tet.opp_tets[2] = TetWalker::new(adj0t.id(), adj0t.edge);
        tet.opp_tets[3] = TetWalker::new(id1, 3);

        let tet = &mut mesh.tets[id1];
        tet.vertices = [v1, v0, v2, v4];
        tet.opp_tets[0] = TetWalker::new(adj1b.id(), adj1b.to_nfe().edge);
        tet.opp_tets[1] = TetWalker::new(adj2b.id(), adj2b.to_pfe().edge);
        tet.opp_tets[2] = TetWalker::new(adj0b.id(), adj0b.edge);
        tet.opp_tets[3] = TetWalker::new(id0, 3);

        mesh.tets.remove(other2.id());

        [TetWalker::new(id0, 10), TetWalker::new(id1, 10)]
    }

    /// Returns the boundary of the region of tets that satisfy some predicate, starting at the current tet.
    /// Also returns all tets in the region.
    /// Each item in the boundary is a tet walker whose current triangle is a triangle of the boundary.
    pub fn boundary_and_enclosed<V, T, F: FnMut(TetId) -> bool>(
        self,
        mesh: &Tets<V, T>,
        mut pred: F,
    ) -> (Vec<TetWalker>, FnvHashSet<TetId>) {
        let mut bound_imm = FnvHashSet::default();
        let mut boundary = vec![];
        let mut enclosed = FnvHashSet::default();
        enclosed.insert(self.id());

        // Initialize intermediate boundary
        bound_imm.insert(BoundaryTri(self));
        bound_imm.insert(BoundaryTri(self.to_opp_edge()));
        bound_imm.insert(BoundaryTri(self.to_twin_edge()));
        bound_imm.insert(BoundaryTri(self.to_perm(Permutation::_3210)));

        // Extend intermediate boundary
        while let Some(tri) = bound_imm.iter().next().copied() {
            bound_imm.remove(&tri);

            let adj = tri.0.to_adj_ae(&mesh);
            if pred(adj.id()) {
                enclosed.insert(adj.id());

                // Local extension
                for walker in [
                    adj.to_opp_edge(),
                    adj.to_twin_edge(),
                    adj.to_perm(Permutation::_3210),
                ]
                .iter()
                {
                    if !bound_imm.remove(&BoundaryTri(walker.to_adj_ae(&mesh))) {
                        bound_imm.insert(BoundaryTri(*walker));
                    }
                }
            } else {
                // Reached edge of boundary
                boundary.push(tri.0);
            }
        }

        (boundary, enclosed)
    }
}

/// A manifold tetrahedralization.
#[derive(Debug)]
pub struct Tets<V, T> {
    vertices: IdMap<VertexId, Vertex<V>>,
    tets: IdMap<TetId, Tet<T>>,
    default_tet: fn() -> T,
}

impl<V, T> Tets<V, T> {
    pub const GHOST: VertexId = VertexId::invalid();

    /// Creates a new tetrahedralization from 4 vertices because it takes
    /// 4 vertices to make a tet. Also takes a default value function for a tet.
    /// Returns the tetrahedralization.
    ///
    /// The vertex ids are VertexId(i) for i in 0..4.
    ///
    /// The solid tet id is TetId(0) and the ghost tet ids are TetId(i) for i in 1..5.
    fn new(vertices: [(Pt3, V); 4], default_tet: fn() -> T) -> Self
    where
        V: Clone,
        T: Clone,
    {
        let mut tets = Self {
            vertices: IdMap::default(),
            tets: IdMap::default(),
            default_tet,
        };

        // Make sure the tet is oriented positive
        let swap = !sim::orient_3d(&vertices, |l, i| l[i].0.coords, 0, 1, 2, 3);

        let tw = |id, edge| TetWalker::new(TetId(id), edge);
        tets.vertices
            .extend_values(vertices.iter().enumerate().map(|(i, (pos, v))| {
                // I'm adding the solid (not ghost) tet first, so TetId(0) is fine here
                Vertex::new(tw(0, 3 - i as IdType), *pos, v.clone())
            }))
            .for_each(|_| {});
        let vi = [
            VertexId(0),
            VertexId(1),
            VertexId(2 + swap as IdType),
            VertexId(3 - swap as IdType),
            Self::GHOST,
        ];

        let vertices_arr = [
            [vi[0], vi[1], vi[2], vi[3]],
            [vi[1], vi[2], vi[3], vi[4]],
            [vi[0], vi[3], vi[2], vi[4]],
            [vi[3], vi[0], vi[1], vi[4]],
            [vi[2], vi[1], vi[0], vi[4]],
        ];

        let opps_arr = [
            [tw(1, 7), tw(2, 7), tw(3, 7), tw(4, 7)],
            [tw(2, 8), tw(3, 5), tw(4, 2), tw(0, 4)],
            [tw(1, 8), tw(4, 5), tw(3, 2), tw(0, 5)],
            [tw(4, 8), tw(1, 5), tw(2, 2), tw(0, 6)],
            [tw(3, 8), tw(2, 5), tw(1, 2), tw(0, 7)],
        ];

        tets.tets
            .extend_values(
                vertices_arr
                    .iter()
                    .zip(opps_arr.iter())
                    .map(|(vertices, opp_tets)| Tet::new(*vertices, *opp_tets, default_tet())),
            )
            .for_each(|_| {});

        tets
    }

    /// Gets the number of vertices in the tet mesh.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets the number of tets in the tet mesh.
    pub fn num_tets(&self) -> usize {
        self.tets.len()
    }

    /// Gets a vertex, if it exists and is solid (not ghost).
    pub fn vertex(&self, vertex: VertexId) -> Option<&Vertex<V>> {
        self.vertices.get(vertex)
    }

    /// Gets a tet, if it exists. It may be a ghost tet.
    pub fn tet(&self, tet: TetId) -> Option<&Tet<T>> {
        self.tets.get(tet)
    }

    /// Gets a tet walker that starts at the given vertex.
    /// The vertex must be a solid vertex.
    pub fn walker_from_vertex(&self, vertex: VertexId) -> TetWalker {
        self.vertices[vertex].tet
    }

    /// Gets a canonicalized tet walker that starts at the given tet.
    /// The tet must exist
    pub fn walker_from_tet(&self, tet: TetId) -> TetWalker {
        TetWalker::new(tet, 3)
    }
}

impl<V, T> Index<VertexId> for Tets<V, T> {
    type Output = Vertex<V>;

    fn index(&self, index: VertexId) -> &Self::Output {
        self.vertex(index)
            .unwrap_or_else(|| panic!("Vertex {} does not exist", index))
    }
}

impl<V, T> Index<TetId> for Tets<V, T> {
    type Output = Tet<T>;

    fn index(&self, index: TetId) -> &Self::Output {
        self.tet(index)
            .unwrap_or_else(|| panic!("Tet {} does not exist", index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! sorted_fn {
        ($name:ident, $n:expr) => {
            fn $name<Idx: Ord + Copy>(mut arr: [Idx; $n]) -> [Idx; $n] {
                arr.sort();
                arr
            }
        };
    }

    macro_rules! even_sort_fn {
        ($name:ident, $n:expr) => {
            /// Even-sorts an array of $n elements
            /// and returns the even-sorted array.
            ///
            /// An even sort sorts the list, but leaves
            /// the last 2 elements out of order if
            /// the original permutation was odd.
            fn $name<Idx: Ord + Copy>(mut arr: [Idx; $n]) -> [Idx; $n] {
                let mut num_swaps = 0;

                for i in 1..$n {
                    for j in (0..i).rev() {
                        if arr[j] > arr[j + 1] {
                            arr.swap(j, j + 1);
                            num_swaps += 1;
                        } else {
                            break;
                        }
                    }
                }

                if num_swaps % 2 != 0 {
                    let len = arr.len();
                    arr.swap(len - 2, len - 1);
                }
                arr
            }
        };
    }

    sorted_fn!(sorted_3, 3);
    even_sort_fn!(even_sort_3, 3);
    even_sort_fn!(even_sort_4, 4);

    fn v(n: IdType) -> VertexId {
        VertexId(n)
    }

    fn t(n: IdType) -> TetId {
        TetId(n)
    }

    fn default_tets() -> Tets<(), ()> {
        Tets::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
            ],
            || (),
        )
    }

    fn reverse_tets() -> Tets<(), ()> {
        Tets::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
            ],
            || (),
        )
    }

    #[test]
    fn test_new() {
        let mesh = default_tets();
        assert_eq!(mesh.num_vertices(), 4);
        assert_eq!(mesh[v(0)].position(), Pt3::new(0.0, 0.0, 0.0));
        assert_eq!(mesh[v(1)].position(), Pt3::new(0.0, 0.0, 1.0));
        assert_eq!(mesh[v(2)].position(), Pt3::new(0.0, 1.0, 0.0));
        assert_eq!(mesh[v(3)].position(), Pt3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.vertex(Tets::<(), ()>::GHOST).map(|v| v.value()), None);

        assert_eq!(mesh.num_tets(), 5);
        assert_eq!(even_sort_4(mesh[t(0)].vertices()), [v(0), v(1), v(2), v(3)]);
        assert!(mesh[t(1)].vertices().contains(&Tets::<(), ()>::GHOST));
    }

    #[test]
    fn test_inverted() {
        // Vertices given in negative orientation order.
        let mesh = reverse_tets();
        assert_eq!(mesh.num_vertices(), 4);
        assert_eq!(mesh[v(0)].position(), Pt3::new(0.0, 0.0, 0.0));
        assert_eq!(mesh[v(1)].position(), Pt3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh[v(2)].position(), Pt3::new(0.0, 1.0, 0.0));
        assert_eq!(mesh[v(3)].position(), Pt3::new(0.0, 0.0, 1.0));
        assert_eq!(mesh.vertex(Tets::<(), ()>::GHOST).map(|v| v.value()), None);

        assert_eq!(mesh.num_tets(), 5);
        assert_eq!(even_sort_4(mesh[t(0)].vertices()), [v(0), v(1), v(3), v(2)]);
        assert!(mesh[t(1)].vertices().contains(&Tets::<(), ()>::GHOST));
    }

    #[test]
    fn test_walker_from_vertex() {
        let mesh = default_tets();

        for i in 0..4 {
            let walker = mesh.walker_from_vertex(v(i));
            assert_eq!(walker.first(&mesh), v(i));
        }
    }

    #[test]
    fn test_walker_canon() {
        let mesh = default_tets();
        let walker = mesh.walker_from_vertex(v(1)).to_canon_tet();
        assert_eq!(walker.tet(&mesh), mesh[walker.id()].vertices());
    }

    #[test]
    fn test_walker() {
        let mesh = default_tets();
        let vs = mesh[t(0)].vertices();
        let walker = mesh.walker_from_tet(t(0));

        assert_eq!(walker.first(&mesh), vs[0]);
        assert_eq!(walker.second(&mesh), vs[1]);
        assert_eq!(walker.third(&mesh), vs[2]);
        assert_eq!(walker.fourth(&mesh), vs[3]);
        assert_eq!(walker.edge(&mesh), [vs[0], vs[1]]);
        assert_eq!(walker.opp_edge(&mesh), [vs[2], vs[3]]);
        assert_eq!(walker.tri(&mesh), [vs[0], vs[1], vs[2]]);
        assert_eq!(walker.opp_tri(&mesh), [vs[3], vs[2], vs[1]]);
        assert_eq!(walker.tet(&mesh), [vs[0], vs[1], vs[2], vs[3]]);
        assert_eq!(walker.id(), t(0));

        // Make sure edges are correct
        #[rustfmt::skip]
        {
            assert_eq!(TetWalker::new(t(0),  0).tet(&mesh), [vs[3], vs[2], vs[1], vs[0]]);
            assert_eq!(TetWalker::new(t(0),  1).tet(&mesh), [vs[2], vs[3], vs[0], vs[1]]);
            assert_eq!(TetWalker::new(t(0),  2).tet(&mesh), [vs[1], vs[0], vs[3], vs[2]]);
            assert_eq!(TetWalker::new(t(0),  3).tet(&mesh), [vs[0], vs[1], vs[2], vs[3]]);
            assert_eq!(TetWalker::new(t(0),  4).tet(&mesh), [vs[2], vs[1], vs[3], vs[0]]);
            assert_eq!(TetWalker::new(t(0),  5).tet(&mesh), [vs[3], vs[0], vs[2], vs[1]]);
            assert_eq!(TetWalker::new(t(0),  6).tet(&mesh), [vs[0], vs[3], vs[1], vs[2]]);
            assert_eq!(TetWalker::new(t(0),  7).tet(&mesh), [vs[1], vs[2], vs[0], vs[3]]);
            assert_eq!(TetWalker::new(t(0),  8).tet(&mesh), [vs[1], vs[3], vs[2], vs[0]]);
            assert_eq!(TetWalker::new(t(0),  9).tet(&mesh), [vs[0], vs[2], vs[3], vs[1]]);
            assert_eq!(TetWalker::new(t(0), 10).tet(&mesh), [vs[3], vs[1], vs[0], vs[2]]);
            assert_eq!(TetWalker::new(t(0), 11).tet(&mesh), [vs[2], vs[0], vs[1], vs[3]]);
        };

        // Make sure operations are correct
        for i in 0..12 {
            let walker = TetWalker::new(t(0), i);
            let mut vertices = walker.tet(&mesh);

            vertices[..3].rotate_left(1);
            assert_eq!(walker.to_next_tri_edge().tet(&mesh), vertices);

            vertices[..3].rotate_left(1);
            assert_eq!(walker.to_prev_tri_edge().tet(&mesh), vertices);

            vertices[..3].rotate_left(1);
            vertices[1..].rotate_left(1);
            assert_eq!(walker.to_next_vertex_edge().tet(&mesh), vertices);

            vertices[1..].rotate_left(1);
            assert_eq!(walker.to_prev_vertex_edge().tet(&mesh), vertices);

            vertices[1..].rotate_left(1);
            vertices.swap(0, 1);
            vertices.swap(2, 3);
            assert_eq!(walker.to_twin_edge().tet(&mesh), vertices);

            vertices.swap(0, 1);
            vertices.swap(2, 3);
            vertices.rotate_left(2);
            assert_eq!(walker.to_opp_edge().tet(&mesh), vertices);
        }
    }

    #[test]
    fn test_to_twin_tri() {
        for mesh in vec![default_tets(), reverse_tets()] {
            for tet in 0..5 {
                for i in 0..12 {
                    let walker = TetWalker::new(t(tet), i);
                    let mut vertices = walker.tet(&mesh);

                    // to_twin tri should flip the current edge and set the fourth vertex to the opposite one
                    vertices.swap(1, 0);
                    vertices[3] = v((Tets::<(), ()>::GHOST.0 as u64 + 6
                        - vertices.iter().map(|v| v.0 as u64).sum::<u64>())
                        as IdType);
                    assert_eq!(
                        walker.to_adj(&mesh).tet(&mesh),
                        vertices,
                        "Flipping {:?}, expected {:?}, got {:?}",
                        walker.tet(&mesh),
                        vertices,
                        walker.to_adj(&mesh).tet(&mesh)
                    );
                }
            }
        }
    }

    #[test]
    fn test_flip14_unchecked() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let walkers = mesh.walker_from_tet(t(1)).flip14_unchecked(&mut mesh, v(4));
        for (i, walker) in walkers.iter().enumerate() {
            assert_eq!(walker.fourth(&mesh), v(4));
        }

        for tet in 0..8 {
            for i in 0..12 {
                let walker = TetWalker::new(t(tet), i);
                let mut vertices = walker.tet(&mesh);

                if !vertices[..3].contains(&v(0)) {
                    vertices.swap(0, 1);

                    if vertices[3] == v(4) {
                        // Internal to external link
                        vertices[3] = v(0);
                    } else if vertices[3] == v(0) {
                        // External to internal link
                        vertices[3] = v(4);
                    } else {
                        // Internal link
                        vertices[3] = v((Tets::<(), ()>::GHOST.0 as u64 + 10
                            - vertices.iter().map(|v| v.0 as u64).sum::<u64>())
                            as IdType);
                    }
                    assert_eq!(
                        walker.to_adj(&mesh).tet(&mesh),
                        vertices,
                        "Flipping {:?}, expected {:?}, got {:?}",
                        walker.tet(&mesh),
                        vertices,
                        walker.to_adj(&mesh).tet(&mesh)
                    );
                }
            }
        }

        for perm in 0..12 {
            print!("Perm {:2}:", perm);
            for i in 0..12 {
                print!("{:3}", TetWalker::PERMUTATION_FWD[perm][i]);
            }
            println!();
        }
        println!();

        for perm in 0..12 {
            print!("Perm {:2}:", perm);
            for i in 0..12 {
                print!("{:3}", TetWalker::PERMUTATION_INV[perm][i]);
            }
            println!();
        }
    }

    #[test]
    fn test_flip23_unchecked() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let walkers = mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4)); // ids 0, 5, 6, 7
        let walkers = walkers[0].flip23_unchecked(&mut mesh); // flips tri [2, 1, 3] away, opp verts 4 and GHOST
        assert_eq!(walkers[0].opp_edge(&mesh), walkers[1].opp_edge(&mesh));
        assert_eq!(walkers[1].opp_edge(&mesh), walkers[2].opp_edge(&mesh));
        assert!(walkers[0].opp_edge(&mesh).contains(&v(4)));
        assert!(walkers[0].opp_edge(&mesh).contains(&Tets::<(), ()>::GHOST));

        for tet in 0..9 {
            for i in 0..12 {
                let walker = TetWalker::new(t(tet), i);
                let mut vertices = walker.tet(&mesh);

                // Resultant tet or opposite resultant tet
                if vertices[3] == v(0)
                    || (vertices.contains(&v(4)) && vertices.contains(&Tets::<(), ()>::GHOST))
                {
                    vertices.swap(0, 1);

                    if vertices[3] == v(4) || vertices[3] == Tets::<(), ()>::GHOST {
                        // Internal to external link
                        vertices[3] = v(0);
                    } else if vertices[3] == v(0) {
                        // External to internal link
                        vertices[3] = if vertices[..3].contains(&v(4)) {
                            Tets::<(), ()>::GHOST
                        } else {
                            v(4)
                        };
                    } else {
                        // Internal link
                        vertices[3] = v((Tets::<(), ()>::GHOST.0 as u64 + 10
                            - vertices.iter().map(|v| v.0 as u64).sum::<u64>())
                            as IdType);
                    }
                    assert_eq!(
                        walker.to_adj(&mesh).tet(&mesh),
                        vertices,
                        "Flipping {:?}, expected {:?}, got {:?}",
                        walker.tet(&mesh),
                        vertices,
                        walker.to_adj(&mesh).tet(&mesh)
                    );
                }
            }
        }
    }

    #[test]
    fn test_flip32_unchecked() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let walkers = mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4)); // ids 0, 5, 6, 7
        let walkers = walkers[0].flip23_unchecked(&mut mesh); // flips tri [2, 1, 3] away, opp verts 4 and GHOST
        let walkers = walkers[0]
            .to_perm(Permutation::_3210)
            .flip32_unchecked(&mut mesh); // Flips in tri [2, 1, 3]
        assert_eq!(sorted_3(walkers[0].opp_tri(&mesh)), [v(1), v(2), v(3)]);
        assert_eq!(sorted_3(walkers[1].opp_tri(&mesh)), [v(1), v(2), v(3)]);
        assert_ne!(
            even_sort_3(walkers[0].opp_tri(&mesh)),
            even_sort_3(walkers[1].opp_tri(&mesh))
        );
        assert_eq!(walkers[0].fourth(&mesh), walkers[1].fourth(&mesh));

        for tet in mesh.tets.keys() {
            for i in 0..12 {
                let walker = TetWalker::new(tet, i);
                let mut vertices = walker.tet(&mesh);

                // Resultant tet or opposite resultant tet
                if vertices[3] == v(0)
                    || (vertices.contains(&v(1))
                        && vertices.contains(&v(2))
                        && vertices.contains(&v(3)))
                {
                    vertices.swap(0, 1);

                    if vertices[3] == v(1) || vertices[3] == v(2) || vertices[3] == v(3) {
                        // Internal to external link
                        vertices[3] = v(0);
                    } else if vertices[3] == v(0) {
                        // External to internal link
                        vertices[3] = v(6 - vertices[..3]
                            .iter()
                            .map(|v| v.0)
                            .filter(|n| *n == 1 || *n == 2 || *n == 3)
                            .sum::<u32>());
                    } else {
                        // Internal link
                        vertices[3] = if vertices[3] == v(4) {
                            Tets::<(), ()>::GHOST
                        } else {
                            v(4)
                        }
                    }
                    assert_eq!(
                        walker.to_adj(&mesh).tet(&mesh),
                        vertices,
                        "Flipping {:?}, expected {:?}, got {:?}",
                        walker.tet(&mesh),
                        vertices,
                        walker.to_adj(&mesh).tet(&mesh)
                    );
                }
            }
        }
    }

    #[test]
    fn test_boundary_and_enclosure_all() {
        let mut mesh = default_tets();
        // Just care about combinatorial structure for this test
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let (boundary, enclosed) = mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4))[0]
            .boundary_and_enclosed(&mesh, |_| true);

        assert_eq!(boundary, vec![]);
        assert_eq!(enclosed, mesh.tets.keys().collect::<FnvHashSet<_>>());
    }

    #[test]
    fn test_boundary_and_enclosure_one() {
        let mut mesh = default_tets();
        // Just care about combinatorial structure for this test
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));

        let walker = mesh.walker_from_tet(t(0));
        let (boundary, enclosed) = walker.flip14_unchecked(&mut mesh, v(4))[0]
            .boundary_and_enclosed(&mesh, |t| mesh[t].vertices() == mesh[walker.id()].vertices());

        assert_eq!(boundary.into_iter().map(BoundaryTri).collect::<FnvHashSet<_>>(), vec![
            BoundaryTri(walker),
            BoundaryTri(walker.to_twin_edge()),
            BoundaryTri(walker.to_opp_edge()),
            BoundaryTri(walker.to_perm(Permutation::_3210)),
        ].into_iter().collect::<FnvHashSet<_>>());

        assert_eq!(enclosed, std::iter::once(walker.id()).collect::<FnvHashSet<_>>());
    }

    #[test]
    fn test_boundary_and_enclosure_some() {
        let mut mesh = default_tets();
        // Just care about combinatorial structure for this test
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));

        let walker = mesh.walker_from_tet(t(0));
        let (boundary, enclosed) = walker.flip14_unchecked(&mut mesh, v(4))[0]
            .boundary_and_enclosed(&mesh, |t| !mesh[t].vertices().contains(&Tets::<(), ()>::GHOST));

        assert_eq!(boundary.into_iter().map(|w| even_sort_3(w.tri(&mesh))).collect::<FnvHashSet<_>>(), vec![
            [v(0), v(1), v(2)],
            [v(1), v(3), v(2)],
            [v(0), v(2), v(3)],
            [v(0), v(3), v(1)],
        ].into_iter().collect::<FnvHashSet<_>>());

        assert_eq!(enclosed, vec![t(0), t(5), t(6), t(7)].into_iter().collect::<FnvHashSet<_>>());
    }
}
