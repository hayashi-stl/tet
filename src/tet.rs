use std::{cell::Cell, collections::VecDeque, hash::{Hash, Hasher}, path::Path};
use std::ops::Index;
use std::iter::{self, Map};
use float_ord::FloatOrd;
use bitflags::bitflags;

use crate::{Pt1, id_map::{self, IdMap, IdType}};
use crate::{Pt3, VertexId, Vec3};
use crate::util;
use fnv::{FnvHashMap, FnvHashSet};
use simplicity as sim;

crate::id! {
    /// A tet mesh tet id
    pub struct TetId
}

bitflags! {
    struct VertexFlags: u32 {
        /// Whether this tet has been visited by vertex_targets_opt
        const VERTEX_TARGETS = 1 << 0;
    }
}

/// A vertex stores a tet that it's part of and its position.
#[derive(Clone, Debug)]
pub struct Vertex<V> {
    tet: TetWalker,
    flags: Cell<VertexFlags>,
    position: Pt3,
    value: V,
}

impl<V> Vertex<V> {
    fn new(tet: TetWalker, position: Pt3, value: V) -> Self {
        Self {
            tet,
            flags: Cell::new(VertexFlags::empty()),
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

bitflags! {
    struct TetFlags: u32 {
        /// Part of the intermediate boundary in boundary_and_enclosed
        const BOUND_IMM_0 = 1 << 0;
        const BOUND_IMM_1 = 1 << 1;
        const BOUND_IMM_2 = 1 << 2;
        const BOUND_IMM_3 = 1 << 3;
        /// Whether this tet is part of the boundary. Used in flip_in_vertex_unchecked.
        const BOUNDARY = 1 << 4;
        /// Whether this tet has been visited by walkers_from_vertex_opt
        const WALKERS_FROM_VERTEX = 1 << 5;
        /// Whether this tet needs to be duplicated in flip_in_vertex_unchecked
        const NEEDS_DUPLICATION = 1 << 6;
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
    flags: Cell<TetFlags>,
    value: T,
}

impl<T> Tet<T> {
    fn new(vertices: [VertexId; 4], opp_tets: [TetWalker; 4], value: T) -> Self {
        Self {
            vertices,
            opp_tets,
            flags: Cell::new(TetFlags::empty()),
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

    fn set_flags(&self, flag: TetFlags) {
        self.flags.set(self.flags.get() | flag);
    }

    fn clear_flags(&self, flag: TetFlags) {
        self.flags.set(self.flags.get() & !flag);
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
    pub fn first<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2][self.edge as usize]]
    }

    /// Gets the second vertex of the tet walker. This is the current edge's target.
    pub fn second<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[2, 3, 0, 1, 1, 0, 3, 2, 3, 2, 1, 0][self.edge as usize]]
    }

    /// Gets the third vertex of the tet walker.
    pub fn third<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1][self.edge as usize]]
    }

    /// Gets the fourth vertex of the tet walker. This is the vertex opposite the current triangle.
    pub fn fourth<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3][self.edge as usize]]
    }

    /// Gets the current edge of the tet walker.
    pub fn edge<V, T>(self, mesh: &TetMesh<V, T>) -> [VertexId; 2] {
        [self.first(mesh), self.second(mesh)]
    }

    /// Gets the opposite edge of the current edge on the tet walker's tet.
    pub fn opp_edge<V, T>(self, mesh: &TetMesh<V, T>) -> [VertexId; 2] {
        [self.third(mesh), self.fourth(mesh)]
    }

    /// Gets the current triangle of the tet walker.
    pub fn tri<V, T>(self, mesh: &TetMesh<V, T>) -> [VertexId; 3] {
        [self.first(mesh), self.second(mesh), self.third(mesh)]
    }

    /// Gets the opposite triangle of the current vertex on the tet walker's tet.
    /// The vertices are returned in [fourth, third, second] order.
    pub fn opp_tri<V, T>(self, mesh: &TetMesh<V, T>) -> [VertexId; 3] {
        [self.fourth(mesh), self.third(mesh), self.second(mesh)]
    }

    /// Gets the vertices of the current tet, respecting the orientation.
    pub fn tet<V, T>(self, mesh: &TetMesh<V, T>) -> [VertexId; 4] {
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
        pub fn to_twin_tri<V, T>((mut) self: Self, mesh: &TetMesh<V, T>) -> Self {
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
        pub fn to_twin_tri_any_edge<V, T>(self: Self, mesh: &TetMesh<V, T>) -> Self {
            mesh[self.id()].opp_tets[self.edge as usize % 4]
        }
    }

    /// Sets the orientation such that the vertices returned by `self.vertices()`
    /// are in the same order as those returned by `mesh[self.id()].vertices()`.
    pub fn to_canon_tet(mut self) -> Self {
        self.edge = 3;
        self
    }

    /// Goes to the canonical orientation of the current triangle.
    pub fn to_canon_tri(mut self) -> Self {
        self.edge %= 4;
        self
    }

    /// Performs a 1-to-4 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Returns walkers of the 4 tets added. Each walker's 4th vertex is the new vertex.
    pub fn flip14_unchecked<V: Clone, T: Clone>(
        self,
        mesh: &mut TetMesh<V, T>,
        vertex: VertexId,
    ) -> [TetWalker; 4] {
        let id0 = self.id();
        let id1 = mesh.tets.insert(mesh.tets[id0].clone());
        let id2 = mesh.tets.insert(mesh.tets[id0].clone());
        let id3 = mesh.tets.insert(mesh.tets[id0].clone());

        // The vertex didn't have a tet before.
        mesh.vertices[vertex].tet = self;

        // Just in case of dangling vertex-tet reference
        let vs = mesh.tets[id0].vertices();
        mesh.vertices.get_mut(vs[0]).map(|v| v.tet = TetWalker::new(id2, 3) );
        mesh.vertices.get_mut(vs[1]).map(|v| v.tet = TetWalker::new(id3, 2) );
        mesh.vertices.get_mut(vs[2]).map(|v| v.tet = TetWalker::new(id0, 1) );
        mesh.vertices.get_mut(vs[3]).map(|v| v.tet = TetWalker::new(id1, 0) );

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
    pub fn flip23_unchecked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> [TetWalker; 3] {
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

        // Just in case of dangling vertex-tet reference
        mesh.vertices.get_mut(v0).map(|v| v.tet = TetWalker::new(id0, 3) );
        mesh.vertices.get_mut(v1).map(|v| v.tet = TetWalker::new(id1, 3) );
        mesh.vertices.get_mut(v2).map(|v| v.tet = TetWalker::new(id2, 3) );
        mesh.vertices.get_mut(v3).map(|v| v.tet = TetWalker::new(id0, 1) );
        mesh.vertices.get_mut(v4).map(|v| v.tet = TetWalker::new(id0, 0) );

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
    pub fn flip32_unchecked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> [TetWalker; 2] {
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

        // Just in case of dangling vertex-tet reference
        mesh.vertices.get_mut(v0).map(|v| v.tet = TetWalker::new(id0, 3) );
        mesh.vertices.get_mut(v1).map(|v| v.tet = TetWalker::new(id0, 7) );
        mesh.vertices.get_mut(v2).map(|v| v.tet = TetWalker::new(id0, 11) );
        mesh.vertices.get_mut(v3).map(|v| v.tet = TetWalker::new(id0, 0) );
        mesh.vertices.get_mut(v4).map(|v| v.tet = TetWalker::new(id1, 0) );

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
        mesh: &TetMesh<V, T>,
        mut pred: F,
    ) -> (Vec<TetWalker>, Vec<TetId>) {
        // Initialize intermediate boundary
        let mut bound_imm = vec![
            TetWalker::new(self.id(), 0),
            TetWalker::new(self.id(), 1),
            TetWalker::new(self.id(), 2),
            TetWalker::new(self.id(), 3),
        ];
        mesh.tets[self.id()].set_flags(TetFlags::BOUND_IMM_0 | TetFlags::BOUND_IMM_1 | TetFlags::BOUND_IMM_2 | TetFlags::BOUND_IMM_3);

        let mut boundary = vec![];
        let mut enclosed = vec![self.id()];

        let edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::BOUND_IMM_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::BOUND_IMM_0.bits << edge % 4)
        };

        // Extend intermediate boundary
        while let Some(tri) = bound_imm.pop() {
            // Clear flag
            let tet = &mesh.tets[tri.id()];
            if !tet.flags.get().intersects(edge_flag(tri.edge)) {
                // Triangle was removed by twin triangle.
                continue;
            }
            tet.clear_flags(edge_flag(tri.edge));

            let adj = tri.to_adj_ae(&mesh);
            if pred(adj.id()) {
                // Each tet should be visited at most once.
                enclosed.push(adj.id());

                // Local extension
                for walker in [
                    adj.to_opp_edge(),
                    adj.to_twin_edge(),
                    adj.to_perm(Permutation::_3210),
                ]
                .iter()
                {
                    let twin = walker.to_adj_ae(&mesh);
                    if mesh.tets[twin.id()].flags.get().intersects(edge_flag(twin.edge)) {
                        mesh.tets[twin.id()].clear_flags(edge_flag(twin.edge));
                    } else {
                        mesh.tets[walker.id()].set_flags(edge_flag(walker.edge));
                        bound_imm.push(*walker);
                    }
                }
            } else {
                // Reached edge of boundary
                boundary.push(tri);
            }
        }

        debug_assert!(
            enclosed.len() == enclosed.iter().collect::<FnvHashSet<_>>().len(),
            "Enclosure contains repeat elements unexpectedly."
        );
        (boundary, enclosed)
    }
}

/// A manifold tetrahedralization.
#[derive(Debug)]
pub struct TetMesh<V, T> {
    vertices: IdMap<VertexId, Vertex<V>>,
    ghost_flags: Cell<VertexFlags>,
    tets: IdMap<TetId, Tet<T>>,
    default_tet: fn() -> T,
}

impl<V, T> TetMesh<V, T> {
    pub const GHOST: VertexId = VertexId::invalid();

    /// Creates a new tetrahedralization from 4 vertices because it takes
    /// 4 vertices to make a tet. Also takes a default value function for a tet.
    /// Returns the tetrahedralization.
    ///
    /// The vertex ids are VertexId(i) for i in 0..4.
    ///
    /// The solid tet id is TetId(0) and the ghost tet ids are TetId(i) for i in 1..5.
    pub fn new(vertices: [(Pt3, V); 4], default_tet: fn() -> T) -> Self
    where
        V: Clone,
        T: Clone,
    {
        Self::with_ids([
            (VertexId(0), vertices[0].0, vertices[0].1.clone()),
            (VertexId(1), vertices[1].0, vertices[1].1.clone()),
            (VertexId(2), vertices[2].0, vertices[2].1.clone()),
            (VertexId(3), vertices[3].0, vertices[3].1.clone()),
        ], default_tet)
    }

    fn with_ids(vertices: [(VertexId, Pt3, V); 4], default_tet: fn() -> T) -> Self
    where
        V: Clone,
        T: Clone,
    {
        let mut tets = Self {
            vertices: IdMap::default(),
            ghost_flags: Cell::new(VertexFlags::empty()),
            tets: IdMap::default(),
            default_tet,
        };

        // Make sure the tet is oriented positive
        let swap = !sim::orient_3d(&vertices, |l, i| l[i].1.coords, 0, 1, 2, 3);

        let tw = |id, edge| TetWalker::new(TetId(id), edge);
        tets.vertices
            .extend(vertices.iter().enumerate().map(|(i, (id, pos, v))| {
                // I'm adding the solid (not ghost) tet first, so TetId(0) is fine here
                (*id, Vertex::new(tw(0, 3 - i as IdType), *pos, v.clone()))
            }));
        let vi = [
            vertices[0].0,
            vertices[1].0,
            vertices[2 + swap as usize].0,
            vertices[3 - swap as usize].0,
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

    /// Gets a reference to the flags of the vertex.
    /// If the ghost vertex is passed, this gets the ghost flags.
    fn vertex_flags(&self, vertex: VertexId) -> &Cell<VertexFlags> {
        if vertex == Self::GHOST {
            &self.ghost_flags
        } else {
            &self[vertex].flags
        }
    }

    /// Gets the number of vertices in the tet mesh.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets the number of tets in the tet mesh.
    pub fn num_tets(&self) -> usize {
        self.tets.len()
    }

    /// Iterates over vertex ids and vertices, not including the ghost vertex.
    pub fn vertices(&self) -> Vertices<V> {
        self.vertices.iter()
    }

    /// Iterates over tet ids and tets, including ghost tets.
    pub fn tets(&self) -> Tets<T> {
        self.tets.iter()
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

    /// Gets the tet walkers that start from a vertex.
    /// Each tet walker is for a unique tet and has the vertex as its first vertex.
    pub fn walkers_from_vertex<'a>(&'a self, vertex: VertexId) -> WalkersFromVertex<'a, V, T> {
        WalkersFromVertex {
            mesh: self,
            visited: FnvHashSet::default(),
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the tet walkers that start from a vertex.
    /// Each tet walker is for a unique tet and has the vertex as its first vertex.
    ///
    /// # Safety
    /// The returned iterator needs to drop before calling this function again.
    unsafe fn walkers_from_vertex_opt<'a>(&'a self, vertex: VertexId) -> WalkersFromVertexOpt<'a, V, T> {
        WalkersFromVertexOpt {
            mesh: self,
            visited: vec![],
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the tets that contain a vertex.
    pub fn vertex_tets<'a>(&'a self, vertex: VertexId) -> VertexTets<'a, V, T> {
        self.walkers_from_vertex(vertex).map(|walker| walker.id())
    }

    /// Gets the tets that contain a vertex.
    ///
    /// # Safety
    /// The returned iterator needs to drop before calling this function again.
    unsafe fn vertex_tets_opt<'a>(&'a self, vertex: VertexId) -> VertexTetsOpt<'a, V, T> {
        self.walkers_from_vertex_opt(vertex).map(|walker| walker.id())
    }

    /// Gets the edge target vertices from a vertex.
    pub fn vertex_targets<'a>(&'a self, vertex: VertexId) -> VertexTargets<'a, V, T> {
        VertexTargets {
            mesh: self,
            visited: FnvHashSet::default(),
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the edge target vertices from a vertex.
    ///
    /// # Safety
    /// The returned iterator needs to drop before calling this function again.
    unsafe fn vertex_targets_opt<'a>(&'a self, vertex: VertexId) -> VertexTargetsOpt<'a, V, T> {
        VertexTargetsOpt {
            mesh: self,
            visited: vec![],
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the position of a vertex. The ghost vertex is at infinity.
    pub fn vertex_position(&self, vertex: VertexId) -> Pt3 {
        if vertex == Self::GHOST {
            Pt1::new(f64::INFINITY).xxx()
        } else {
            self[vertex].position()
        }
    }

    /// Gets the distance between two vertices.
    pub fn vertex_distance(&self, v0: VertexId, v1: VertexId) -> f64 {
        self.vertex_distance2(v0, v1).sqrt()
    }

    /// Gets the distance squared between two vertices.
    /// If exactly 1 vertex is the ghost vertex, this returns infinity.
    /// If both vertices are the ghost vertex, this returns 0.
    pub fn vertex_distance2(&self, v0: VertexId, v1: VertexId) -> f64 {
        if v0 == Self::GHOST {
            if v1 == Self::GHOST { 0.0 } else { f64::INFINITY }
        } else if v1 == Self::GHOST {
            f64::INFINITY
        } else {
            (self[v0].position() - self[v1].position()).norm_squared()
        }
    }

    /// Gets the tets adjacent to this tet.
    pub fn adjacent_tets(&self, tet: TetId) -> [TetId; 4] {
        [
            self[tet].opp_tets[0].id(),
            self[tet].opp_tets[1].id(),
            self[tet].opp_tets[2].id(),
            self[tet].opp_tets[3].id(),
        ]
    }

    fn idx(&self, vertex: VertexId) -> Vec3 {
        self[vertex].position().coords
    }

    /// Gets whether the last point is in the circumsphere of the first 4 points.
    /// Uses simulation of simplicity to avoid ties.
    /// [v0, v1, v2, v3] must have positive orientation.
    /// Assumes that none of the points equal and that the last point isn't the ghost vertex.
    pub fn in_sphere(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId, v4: VertexId) -> bool {
        if v0 == Self::GHOST {
            sim::orient_3d(self, Self::idx, v3, v2, v1, v4)
        } else if v1 == Self::GHOST {
            sim::orient_3d(self, Self::idx, v2, v3, v0, v4)
        } else if v2 == Self::GHOST {
            sim::orient_3d(self, Self::idx, v1, v0, v3, v4)
        } else if v3 == Self::GHOST {
            sim::orient_3d(self, Self::idx, v0, v1, v2, v4)
        } else {
            sim::in_sphere(self, Self::idx, v0, v1, v2, v3, v4)
        }
    }

    /// Flips in a vertex using one big mega flip.
    /// The boundary must be manifold.
    /// For now, no vertex can be completely enclosed.
    /// This does not check that negative-volume tets won't be created. Be careful.
    pub fn flip_in_vertex_unchecked(&mut self, vertex: VertexId, boundary: &mut [TetWalker], enclosure: &[TetId]) where T: Clone {
        // Delete tets that don't even have a triangle on the boundary
        // Mark boundary
        for walker in &*boundary {
            self.tets[walker.id()].set_flags(TetFlags::BOUNDARY);
        }
        // Removed unmarked enclosure
        for id in enclosure {
            if !self.tets[*id].flags.get().intersects(TetFlags::BOUNDARY) {
                self.tets.remove(*id);
            }
        }

        // Check for tets that need to be duplicated and make necessary clones
        for walker in boundary.iter_mut() {
            if self.tets[walker.id()].flags.get().intersects(TetFlags::NEEDS_DUPLICATION) {
                let id = self.tets.insert(self[walker.id()].clone());
                // Do not clear the flag in the clone because it will be cleared later.
                let adj = walker.to_adj_ae(&self);
                self.tets[adj.id()].opp_tets[adj.edge as usize % 4].tet = id;
                walker.tet = id;
            } else {
                self.tets[walker.id()].set_flags(TetFlags::NEEDS_DUPLICATION);
            }
        }

        // Build adjacency map; will need for internal links.
        // This is why the boundary must be manifold.
        let mut edge_map = FnvHashMap::<[VertexId; 2], TetWalker>::default();
        edge_map.reserve(3 * boundary.len());

        for walker in &*boundary {
            // Also unmark boundary and tet duplication flags
            self.tets[walker.id()].clear_flags(TetFlags::BOUNDARY | TetFlags::NEEDS_DUPLICATION);

            // Don't forget the new vertex!
            self.tets[walker.id()].vertices[walker.edge as usize % 4] = vertex;
            if self[vertex].tet.id() == TetId::invalid() {
                self.vertices[vertex].tet = walker.to_perm(Permutation::_3210);
            }
            // Other vertices, just in case their references are now dangling
            let tri = walker.tri(self);
            self.vertices.get_mut(tri[0]).map(|v| v.tet = *walker);
            self.vertices.get_mut(tri[1]).map(|v| v.tet = walker.to_nfe());
            self.vertices.get_mut(tri[2]).map(|v| v.tet = walker.to_pfe());

            let tri = walker.tri(&self);
            edge_map.insert([tri[0], tri[1]], *walker);
            edge_map.insert([tri[1], tri[2]], walker.to_nfe());
            edge_map.insert([tri[2], tri[0]], walker.to_pfe());
        }

        // Internal-external links didn't need changing.
        // External-internal links got fixed.
        // Now fix internal links.
        for ([v0, v1], walker) in &edge_map {
            let walker = walker.to_twin_edge();
            let adj = edge_map[&[*v1, *v0]].to_twin_edge();

            // Calling `walker.to_adj(&self)` should give `adj`
            self.tets[walker.id()].opp_tets[walker.edge as usize % 4] = match walker.edge / 4 {
                0 => adj,
                1 => adj.to_nfe(),
                2 => adj.to_pfe(),
                _ => unreachable!(),
            };
        }
    }

    /// Create a Delaunay tetrahedralization from vertices. Panics if there are less than 4 vertices
    /// because it takes 4 vertices to make a tet.
    /// `VertexId(i)` is the index to the `i`th vertex returned from the iterator.
    pub fn delaunay_from_vertices<I: IntoIterator<Item = (Pt3, V)>>(vertices: I, default_tet: fn() -> T) -> Self where V: Clone, T: Clone {
        let mut vertices = util::hilbert_sort(vertices);
        let mut drain = vertices.drain(0..4);
        let v0 = drain.next().unwrap();
        let v1 = drain.next().unwrap();
        let v2 = drain.next().unwrap();
        let v3 = drain.next().unwrap();
        std::mem::drop(drain);
        let v3_id = v3.0;

        let mut mesh = Self::with_ids([v0, v1, v2, v3], default_tet);
        // include v3's id to avoid indexing out of bounds
        let ids = iter::once(v3_id).chain(vertices.iter().map(|v| v.0)).collect::<Vec<_>>();

        for (i, (id, position, value)) in vertices.into_iter().enumerate() {
            mesh.add_vertex_delaunay_internal(Some(id), position, value, ids[i]);
        }

        mesh
    }

    /// Add a vertex to keep the Delaunay property, assuming this is a Delaunay tetrahedralization.
    /// Returns the new vertex id.
    pub fn add_vertex_delaunay(&mut self, position: Pt3, value: V) -> VertexId where T: Clone {
        self.add_vertex_delaunay_internal(None, position, value, self.vertices.keys().next().unwrap())
    }

    pub fn add_vertex_delaunay_internal(&mut self, id: Option<VertexId>, position: Pt3, value: V, search_start: VertexId) -> VertexId where T: Clone {
        let mut closest = search_start;

        let id = if let Some(id) = id {
            self.vertices.insert_with_key(id, Vertex::new(TetWalker::new(TetId::invalid(), 0), position, value));
            id
        } else {
            self.vertices.insert(Vertex::new(TetWalker::new(TetId::invalid(), 0), position, value))
        };

        let lucky_tet = self[closest].tet.id();
        let not_delaunay = if {
            let vs = self[lucky_tet].vertices;
            self.in_sphere(vs[0], vs[1], vs[2], vs[3], id)
        } {
            lucky_tet
        } else {
            // Go to vertex closest to position
            while let Some(closer) = {
                // Safe because this function's iterator is dropped on every iteration.
                let min = unsafe {
                    self.vertex_targets_opt(closest)
                        .filter(|v| *v != id)
                        .min_by_key(|v| FloatOrd(self.vertex_distance2(*v, id)))
                };
                min.filter(|v| self.vertex_distance2(*v, id) < self.vertex_distance2(closest, id))
            } {
                closest = closer;
            }

            // Find a tet that is no longer Delaunay
            // Safe because there was no call to vertex_tets_opt before
            unsafe {
                self.vertex_tets_opt(closest).find(|tet| { 
                    let vs = self[*tet].vertices();
                    self.in_sphere(vs[0], vs[1], vs[2], vs[3], id)
                })
            }.unwrap_or_else(|| {
                // Rare case because the distance comparisons are not robust
                // Safe because the previous call's iterator was dropped.
                let mut checked = unsafe {
                    self.vertex_tets_opt(closest).collect::<FnvHashSet<_>>()
                }; // already checked
                let mut tets = checked.iter().flat_map(|tet| self.adjacent_tets(*tet).to_vec())
                    .collect::<FnvHashSet<_>>().into_iter().collect::<VecDeque<_>>();

                while let Some(tet) = tets.pop_front() {
                    if checked.insert(tet) {
                        let vs = self[tet].vertices();
                        if self.in_sphere(vs[0], vs[1], vs[2], vs[3], id) {
                            return tet;
                        }
                    }
                }

                panic!("Expected some tet to stop being Delaunay when inserting vertex {} at {:?}", id, self[id].position())
            })
        };

        // Find all non-Delaunay tets
        let (mut boundary, enclosed) = self.walker_from_tet(not_delaunay).boundary_and_enclosed(self,
            |tet| {
                let vs = self[tet].vertices();
                self.in_sphere(vs[0], vs[1], vs[2], vs[3], id)
            });

        // Add the new vertex
        self.flip_in_vertex_unchecked(id, &mut boundary, &enclosed);

        id
    }

    #[cfg(feature = "obj")]
    pub fn export_debug_obj<P: AsRef<Path>>(&self, path: P) {
        let solid_tets = || self.tets.values().filter(|t| !t.vertices().contains(&Self::GHOST)).enumerate();

        let obj = obj::ObjData {
            position: solid_tets().flat_map(|(_, t)| t.vertices().to_vec())
                .map(|v| {
                    let pos = self[v].position;
                    [pos.x as f32, pos.y as f32, pos.z as f32]
                }).collect(),
            texture: vec![],
            normal: vec![],

            objects: vec![obj::Object {
                name: "Tet Mesh".to_owned(),
                groups: vec![obj::Group {
                    name: "Tet Group".to_owned(),
                    index: 0,
                    material: None,
                    polys: solid_tets().flat_map(|(i, _)| vec![
                        [4 * i + 0, 4 * i + 1, 4 * i + 2],
                        [4 * i + 3, 4 * i + 2, 4 * i + 1],
                        [4 * i + 2, 4 * i + 3, 4 * i + 0],
                        [4 * i + 1, 4 * i + 0, 4 * i + 3],
                    ]).map(|[a, b, c]| obj::SimplePolygon(vec![
                        obj::IndexTuple(a, None, None),
                        obj::IndexTuple(b, None, None),
                        obj::IndexTuple(c, None, None),
                    ])).collect::<Vec<_>>()
                }],
            }],

            material_libs: vec![],
        };
        obj.save(path).unwrap();
    }
}

impl<V, T> Index<VertexId> for TetMesh<V, T> {
    type Output = Vertex<V>;

    fn index(&self, index: VertexId) -> &Self::Output {
        self.vertex(index)
            .unwrap_or_else(|| panic!("Vertex {} does not exist", index))
    }
}

impl<V, T> Index<TetId> for TetMesh<V, T> {
    type Output = Tet<T>;

    fn index(&self, index: TetId) -> &Self::Output {
        self.tet(index)
            .unwrap_or_else(|| panic!("Tet {} does not exist", index))
    }
}

/// Iterates over vertex ids and vertices.
pub type Vertices<'a, V> = id_map::Iter<'a, VertexId, Vertex<V>>;

/// Iterates over tet ids and tets.
pub type Tets<'a, T> = id_map::Iter<'a, TetId, Tet<T>>;

/// Iterates over tet walkers from a vertex.
/// Each tet walker is for a unique tetrahedron and has the vertex as its first vertex.
#[derive(Clone, Debug)]
pub struct WalkersFromVertex<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: FnvHashSet<TetId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for WalkersFromVertex<'a, V, T> {
    type Item = TetWalker;

    fn next(&mut self) -> Option<Self::Item> {
        // DFS
        while let Some(walker) = self.to_search.pop() {
            if self.visited.insert(walker.id()) {
                self.to_search.push(walker.to_perm(Permutation::_1032).to_adj(self.mesh));
                self.to_search.push(walker.to_perm(Permutation::_2013).to_adj(self.mesh));
                self.to_search.push(walker.to_perm(Permutation::_3021).to_adj(self.mesh));
                return Some(walker)
            }
        }
        None
    }
}

/// Iterates over tet walkers from a vertex.
/// Each tet walker is for a unique tetrahedron and has the vertex as its first vertex.
/// Uses flags and a vec instead of a hash map. Do not call walkers_from_vertex_opt again until this this dropped.
#[derive(Clone, Debug)]
struct WalkersFromVertexOpt<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: Vec<TetId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for WalkersFromVertexOpt<'a, V, T> {
    type Item = TetWalker;

    fn next(&mut self) -> Option<Self::Item> {
        // DFS
        while let Some(walker) = self.to_search.pop() {
            let tet = &self.mesh.tets[walker.id()];
            if !tet.flags.get().intersects(TetFlags::WALKERS_FROM_VERTEX) {
                tet.set_flags(TetFlags::WALKERS_FROM_VERTEX);
                self.visited.push(walker.id());
                self.to_search.push(walker.to_perm(Permutation::_1032).to_adj(self.mesh));
                self.to_search.push(walker.to_perm(Permutation::_2013).to_adj(self.mesh));
                self.to_search.push(walker.to_perm(Permutation::_3021).to_adj(self.mesh));
                return Some(walker)
            }
        }
        None
    }
}

impl<'a, V, T> Drop for WalkersFromVertexOpt<'a, V, T> {
    fn drop(&mut self) {
        for tet in &self.visited {
            self.mesh.tets[*tet].clear_flags(TetFlags::WALKERS_FROM_VERTEX);
        }
    }
}

/// Iterates over tets containing a vertex.
pub type VertexTets<'a, V, T> = Map<WalkersFromVertex<'a, V, T>, fn(TetWalker) -> TetId>;

/// Iterates over tets containing a vertex.
type VertexTetsOpt<'a, V, T> = Map<WalkersFromVertexOpt<'a, V, T>, fn(TetWalker) -> TetId>;

#[derive(Clone, Debug)]
pub struct VertexTargets<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: FnvHashSet<VertexId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for VertexTargets<'a, V, T> {
    type Item = VertexId;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(walker) = self.to_search.pop() {
            let target = walker.second(self.mesh);

            if self.visited.insert(target) {
                // Search is different; search 1 edge radius at a time.
                let mut walker = walker.to_adj(self.mesh).to_nfe();
                let start = walker;
                while {
                    self.to_search.push(walker);
                    walker = walker.to_perm(Permutation::_3021).to_adj(self.mesh);
                    walker != start
                } {}

                return Some(target);
            }
        }
        None
    }
}

#[derive(Clone, Debug)]
pub struct VertexTargetsOpt<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: Vec<VertexId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for VertexTargetsOpt<'a, V, T> {
    type Item = VertexId;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(walker) = self.to_search.pop() {
            let target = walker.second(self.mesh);

            let flags = self.mesh.vertex_flags(target);
            if !flags.get().intersects(VertexFlags::VERTEX_TARGETS) {
                flags.set(flags.get() | VertexFlags::VERTEX_TARGETS);
                self.visited.push(target);

                // Search is different; search 1 edge radius at a time.
                let mut walker = walker.to_adj(self.mesh).to_nfe();
                let start = walker;
                while {
                    self.to_search.push(walker);
                    walker = walker.to_perm(Permutation::_3021).to_adj(self.mesh);
                    walker != start
                } {}

                return Some(target);
            }
        }
        None
    }
}

impl<'a, V, T> Drop for VertexTargetsOpt<'a, V, T> {
    fn drop(&mut self) {
        for target in &self.visited {
            let flags = self.mesh.vertex_flags(*target);
            flags.set(flags.get() & !VertexFlags::VERTEX_TARGETS);
        }
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

    fn default_tets() -> TetMesh<(), ()> {
        TetMesh::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
            ],
            || (),
        )
    }

    fn reverse_tets() -> TetMesh<(), ()> {
        TetMesh::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
            ],
            || (),
        )
    }

    #[track_caller]
    fn assert_tets<V, T, I: IntoIterator<Item = [u32; 4]>>(mesh: &TetMesh<V, T>, tets: I) {
        let result = mesh.tets.iter().map(|(_, tet)| even_sort_4(tet.vertices())).collect::<FnvHashSet<_>>();
        let expect = tets.into_iter().collect::<Vec<_>>();
        assert_eq!(result.len(), expect.len());
        let expect = expect.into_iter().map(|[v0, v1, v2, v3]|
            even_sort_4([v(v0), v(v1), v(v2), v(v3)])).collect::<FnvHashSet<_>>();
        assert_eq!(result, expect);
    }

    #[test]
    fn test_new() {
        let mesh = default_tets();
        assert_eq!(mesh.num_vertices(), 4);
        assert_eq!(mesh[v(0)].position(), Pt3::new(0.0, 0.0, 0.0));
        assert_eq!(mesh[v(1)].position(), Pt3::new(0.0, 0.0, 1.0));
        assert_eq!(mesh[v(2)].position(), Pt3::new(0.0, 1.0, 0.0));
        assert_eq!(mesh[v(3)].position(), Pt3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.vertex(TetMesh::<(), ()>::GHOST).map(|v| v.value()), None);

        assert_eq!(mesh.num_tets(), 5);
        assert_eq!(even_sort_4(mesh[t(0)].vertices()), [v(0), v(1), v(2), v(3)]);
        assert!(mesh[t(1)].vertices().contains(&TetMesh::<(), ()>::GHOST));
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
        assert_eq!(mesh.vertex(TetMesh::<(), ()>::GHOST).map(|v| v.value()), None);

        assert_eq!(mesh.num_tets(), 5);
        assert_eq!(even_sort_4(mesh[t(0)].vertices()), [v(0), v(1), v(3), v(2)]);
        assert!(mesh[t(1)].vertices().contains(&TetMesh::<(), ()>::GHOST));
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
                    vertices[3] = v((TetMesh::<(), ()>::GHOST.0 as u64 + 6
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
        for (_, walker) in walkers.iter().enumerate() {
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
                        vertices[3] = v((TetMesh::<(), ()>::GHOST.0 as u64 + 10
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
        assert!(walkers[0].opp_edge(&mesh).contains(&TetMesh::<(), ()>::GHOST));

        for tet in 0..9 {
            for i in 0..12 {
                let walker = TetWalker::new(t(tet), i);
                let mut vertices = walker.tet(&mesh);

                // Resultant tet or opposite resultant tet
                if vertices[3] == v(0)
                    || (vertices.contains(&v(4)) && vertices.contains(&TetMesh::<(), ()>::GHOST))
                {
                    vertices.swap(0, 1);

                    if vertices[3] == v(4) || vertices[3] == TetMesh::<(), ()>::GHOST {
                        // Internal to external link
                        vertices[3] = v(0);
                    } else if vertices[3] == v(0) {
                        // External to internal link
                        vertices[3] = if vertices[..3].contains(&v(4)) {
                            TetMesh::<(), ()>::GHOST
                        } else {
                            v(4)
                        };
                    } else {
                        // Internal link
                        vertices[3] = v((TetMesh::<(), ()>::GHOST.0 as u64 + 10
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
                            TetMesh::<(), ()>::GHOST
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
        assert_eq!(enclosed.into_iter().collect::<FnvHashSet<_>>(), mesh.tets.keys().collect::<FnvHashSet<_>>());
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

        assert_eq!(enclosed, vec![walker.id()]);
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
            .boundary_and_enclosed(&mesh, |t| !mesh[t].vertices().contains(&TetMesh::<(), ()>::GHOST));

        assert_eq!(boundary.into_iter().map(|w| even_sort_3(w.tri(&mesh))).collect::<FnvHashSet<_>>(), vec![
            [v(0), v(1), v(2)],
            [v(1), v(3), v(2)],
            [v(0), v(2), v(3)],
            [v(0), v(3), v(1)],
        ].into_iter().collect::<FnvHashSet<_>>());

        assert_eq!(enclosed.into_iter().collect::<FnvHashSet<_>>(), vec![t(0), t(5), t(6), t(7)].into_iter().collect::<FnvHashSet<_>>());
    }

    #[test]
    fn test_flip_in_vertex_unchecked_1_4() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let (mut boundary, enclosed) = mesh.walker_from_tet(t(1))
            .boundary_and_enclosed(&mesh, |id| id == t(1));
        mesh.flip_in_vertex_unchecked(v(4), &mut boundary, &enclosed);

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
                        vertices[3] = v((TetMesh::<(), ()>::GHOST.0 as u64 + 10
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
    fn test_flip_in_vertex_unchecked_2_6() {
        let mut mesh = default_tets();
        for _ in 0..2 {
            mesh.vertices.insert(Vertex::new(
                TetWalker::new(TetId::invalid(), 0),
                Pt3::origin(),
                (),
            ));
        }

        mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4));
        let (mut boundary, enclosed) = mesh.walker_from_tet(t(1))
            .boundary_and_enclosed(&mesh, |t| {
                let vs = mesh[t].vertices();
                vs.contains(&v(1)) && vs.contains(&v(2)) && vs.contains(&v(3))
            });

        // 2-to-6 flip.
        mesh.flip_in_vertex_unchecked(v(5), &mut boundary, &enclosed);

        for tet in 0..12 {
            for i in 0..12 {
                let walker = TetWalker::new(t(tet), i);
                let mut vertices = walker.tet(&mesh);

                // Resultant tet or opposite resultant tet
                if vertices[3] == v(0) || vertices.contains(&v(5)) {
                    vertices.swap(0, 1);

                    if vertices[3] == v(5) {
                        // Internal to external link
                        vertices[3] = v(0);
                    } else if vertices[3] == v(0) {
                        // External to internal link
                        vertices[3] = v(5);
                    } else {
                        // Internal link
                        vertices[3] = if vertices[3] == v(4) {
                            TetMesh::<(), ()>::GHOST
                        } else if vertices[3] == TetMesh::<(), ()>::GHOST {
                            v(4)
                        } else {
                            v(6 - vertices
                                .iter()
                                .map(|v| v.0)
                                .filter(|n| *n == 1 || *n == 2 || *n == 3)
                                .sum::<u32>())
                        };
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
    fn test_walkers_from_vertex() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4));

        let walkers = mesh.walkers_from_vertex(v(4)).collect::<Vec<_>>();
        assert_eq!(walkers.len(), 4);
        for walker in walkers {
            assert_eq!(walker.first(&mesh), v(4));
        }
    }

    #[test]
    fn test_vertex_tets() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4));

        let tets = mesh.vertex_tets(v(4)).collect::<Vec<_>>();
        assert_eq!(tets.len(), 4);
        assert_eq!(tets.into_iter().map(|tet| even_sort_4(mesh[tet].vertices())).collect::<FnvHashSet<_>>(), vec![
            [v(1), v(2), v(4), v(3)],
            [v(0), v(2), v(3), v(4)],
            [v(0), v(1), v(4), v(3)],
            [v(0), v(1), v(2), v(4)],
        ].into_iter().collect::<FnvHashSet<_>>());
    }

    #[test]
    fn test_vertex_targets() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4));

        let targets = mesh.vertex_targets(v(4)).collect::<Vec<_>>();
        assert_eq!(targets.len(), 4);
        assert_eq!(targets.into_iter().collect::<FnvHashSet<_>>(), vec![v(0), v(1), v(2), v(3)].into_iter().collect::<FnvHashSet<_>>());
    }

    #[test]
    fn test_in_sphere_ghost() {
        let mut mesh = reverse_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::new(0.5, 0.3, 0.6),
            (),
        ));

        assert!(mesh.in_sphere(TetMesh::<(), ()>::GHOST, v(1), v(2), v(3), v(4)));
        assert!(!mesh.in_sphere(v(0), TetMesh::<(), ()>::GHOST, v(2), v(3), v(4)));
        assert!(!mesh.in_sphere(v(0), v(1), TetMesh::<(), ()>::GHOST, v(3), v(4)));
        assert!(!mesh.in_sphere(v(0), v(1), v(2), TetMesh::<(), ()>::GHOST, v(4)));
    }

    #[test]
    fn test_in_sphere_no_ghost() {
        let mut mesh = reverse_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::new(0.5, 0.3, 0.6),
            (),
        ));

        assert!(mesh.in_sphere(v(0), v(1), v(3), v(2), v(4)));
    }

    #[test]
    fn test_delaunay_tets_single() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(vec![
            (Pt3::new(0.0, 0.0, 0.0), ()),
            (Pt3::new(1.0, 0.0, 0.0), ()),
            (Pt3::new(0.0, 1.0, 0.0), ()),
            (Pt3::new(0.0, 0.0, 1.0), ()),
        ], || ());
        
        assert_eq!(mesh.num_vertices(), 4);
        assert_tets(&mesh, vec![
            [0, 1, 3, 2],
            [3, 1, 0, u32::MAX],
            [1, 3, 2, u32::MAX],
            [0, 2, 3, u32::MAX],
            [2, 0, 1, u32::MAX],
        ]);
    }

    #[test]
    fn test_delaunay_tets_multiple() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(vec![
            (Pt3::new(0.0, 0.0, 0.0), ()),
            (Pt3::new(1.0, 0.0, 0.0), ()),
            (Pt3::new(0.0, 1.0, 0.0), ()),
            (Pt3::new(0.0, 0.0, 1.0), ()),
            (Pt3::new(1.5, 1.5, 1.0), ()),
            (Pt3::new(0.5, 0.5, 0.5), ()),
        ], || ());
        
        assert_eq!(mesh.num_vertices(), 6);
        assert_tets(&mesh, vec![
            [0, 3, 2, 5],
            [0, 1, 3, 5],
            [1, 0, 2, 5],
            [1, 4, 3, 5],
            [3, 4, 2, 5],
            [2, 4, 1, 5],
            [0, 2, 3, u32::MAX],
            [0, 3, 1, u32::MAX],
            [1, 2, 0, u32::MAX],
            [1, 3, 4, u32::MAX],
            [3, 2, 4, u32::MAX],
            [2, 1, 4, u32::MAX],
        ]);
    }

    #[test]
    fn test_delaunay_tets_same_position() {
        // Simulation of simplicity is used. This should be perfectly fine.
        let _mesh = TetMesh::<(), ()>::delaunay_from_vertices(vec![
            (Pt3::new(0.0, 0.0, 0.0), ()),
            (Pt3::new(0.0, 0.0, 0.0), ()),
            (Pt3::new(1.0, 0.0, 0.0), ()),
            (Pt3::new(0.0, 1.0, 0.0), ()),
            (Pt3::new(0.0, 0.0, 1.0), ()),
        ], || ());
    }
}
