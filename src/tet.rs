use bitflags::bitflags;
use float_ord::FloatOrd;
use std::ops::Index;
use std::{
    cell::Cell,
    collections::VecDeque,
    hash::{Hash, Hasher},
    path::Path,
};
use std::{
    hint::unreachable_unchecked,
    iter::{self, Map},
};
use util::CircularListIter;

use crate::{
    id_map::{self, IdMap, IdType},
    Pt1,
};
use crate::{util, Plc};
use crate::{Pt3, Vec3, VertexId};
use fnv::{FnvHashMap, FnvHashSet};
use simplicity as sim;

crate::id! {
    /// An undirected edge id
    struct UndirEdgeId
}

crate::id! {
    /// A tet mesh tet id
    pub struct TetId
}

#[inline(always)]
fn sorted_2(mut arr: [VertexId; 2]) -> [VertexId; 2] {
    if arr[0] > arr[1] {
        arr.swap(0, 1);
    }
    arr
}

fn sorted_3(mut arr: [VertexId; 3]) -> [VertexId; 3] {
    if arr[0] > arr[1] {
        arr.swap(0, 1);
    }
    if arr[1] > arr[2] {
        arr.swap(1, 2);
    }
    if arr[0] > arr[1] {
        arr.swap(0, 1);
    }
    arr
}

fn even_sorted_3(mut arr: [VertexId; 3]) -> [VertexId; 3] {
    if arr[1] < arr[0] && arr[1] < arr[2] {
        arr.swap(0, 1);
        arr.swap(1, 2);
    } else if arr[2] < arr[0] && arr[2] < arr[1] {
        arr.swap(1, 2);
        arr.swap(0, 1);
    }
    arr
}

fn even_sorted_4(mut arr: [VertexId; 4]) -> [VertexId; 4] {
    if arr[3] < arr[0] && arr[3] < arr[1] && arr[3] < arr[2] {
        arr.swap(0, 3);
        arr.swap(1, 2);
        let sub = [arr[1], arr[2], arr[3]];
        arr[1..].copy_from_slice(&even_sorted_3(sub));
    } else {
        let sub = [arr[0], arr[1], arr[2]];
        arr[..3].copy_from_slice(&even_sorted_3(sub));
        let sub = [arr[1], arr[2], arr[3]];
        arr[1..].copy_from_slice(&even_sorted_3(sub));
    }
    arr
}

/// Even-sorted list of 3 vertices
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
struct TriVertices([VertexId; 3]);

impl TriVertices {
    fn new(vertices: [VertexId; 3]) -> Self {
        Self(even_sorted_3(vertices))
    }
}

bitflags! {
    struct VertexFlags: u32 {
        /// Whether this tet has been visited by vertex_targets_opt
        const VERTEX_TARGETS = 1 << 0;
        /// Whether an edge from this vertex to the new vertex has been added.
        /// This is used in flip_in_vertex_unchecked.
        const EDGE_ADDED = 1 << 1;
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

    // There is no set_flags or get_flags method because of the ghost vertex.
}

bitflags! {
    struct UndirEdgeFlags: u32 {
        /// Whether this edge is part of the boundary. Used in flip_in_vertex_unchecked.
        const BOUNDARY = 1 << 0;
        /// Whether this edge is unflippable because it is an edge in the input PLC.
        const LOCKED = 1 << 1;
    }
}

/// An undirected edge just stores its vertices. The vertices are sorted.
#[derive(Clone, Debug)]
struct UndirEdge {
    vertices: [VertexId; 2],
    flags: Cell<UndirEdgeFlags>,
    /// Used by flip_in_vertex_unchecked.
    boundary_walker: TetWalker,
}

impl UndirEdge {
    fn new(vertices: [VertexId; 2]) -> Self {
        Self {
            vertices: sorted_2(vertices),
            flags: Cell::new(UndirEdgeFlags::empty()),
            boundary_walker: TetWalker::new(TetId::invalid(), 0),
        }
    }

    fn set_flags(&self, flag: UndirEdgeFlags) {
        self.flags.set(self.flags.get() | flag);
    }

    fn clear_flags(&self, flag: UndirEdgeFlags) {
        self.flags.set(self.flags.get() & !flag);
    }
}

bitflags! {
    struct TetFlags: u32 {
        /// Part of the intermediate boundary in point_cavity
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
        /// Whether this tet has been removed from the enclosure in point_cavity
        const NOT_ENCLOSURE = 1 << 7;
        /// Whether this tri is unflippable because it is an edge in the input PLC.
        const LOCKED_0 = 1 << 8;
        const LOCKED_1 = 1 << 9;
        const LOCKED_2 = 1 << 10;
        const LOCKED_3 = 1 << 11;
        /// Part of the concave boundary extent in point_cavity
        const BOUND_LIMIT_0 = 1 << 12;
        const BOUND_LIMIT_1 = 1 << 13;
        const BOUND_LIMIT_2 = 1 << 14;
        const BOUND_LIMIT_3 = 1 << 15;
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
    edges: [UndirEdgeId; 6],
    flags: Cell<TetFlags>,
    value: T,
}

impl<T> Tet<T> {
    fn new(
        vertices: [VertexId; 4],
        opp_tets: [TetWalker; 4],
        edges: [UndirEdgeId; 6],
        value: T,
    ) -> Self {
        Self {
            vertices,
            opp_tets,
            flags: Cell::new(TetFlags::empty()),
            edges,
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

/// An edge or triangle
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum EdgeOrTri {
    Edge([VertexId; 2]),
    Tri([VertexId; 3]),
}

/// The type of sliver, either a triangle and a vertex or 2 edges are too close to each other.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SliverKind {
    /// Tet walker on an offending edge.
    EdgeEdge(TetWalker),
    /// Tet walker on the offending triangle.
    TriVertex(TetWalker),
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
    /// Important invariant: this is a valid tet.
    tet: TetId,
    /// Important invariant: `edge` is in 0..12
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
        unsafe {
            mesh.tets.get_unchecked(self.tet).vertices
                [*[3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2].get_unchecked(self.edge as usize)]
        }
    }

    /// Gets the second vertex of the tet walker. This is the current edge's target.
    pub fn second<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        unsafe {
            mesh.tets.get_unchecked(self.tet).vertices
                [*[2, 3, 0, 1, 1, 0, 3, 2, 3, 2, 1, 0].get_unchecked(self.edge as usize)]
        }
    }

    /// Gets the third vertex of the tet walker.
    pub fn third<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        unsafe {
            mesh.tets.get_unchecked(self.tet).vertices
                [*[1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1].get_unchecked(self.edge as usize)]
        }
    }

    /// Gets the fourth vertex of the tet walker. This is the vertex opposite the current triangle.
    pub fn fourth<V, T>(self, mesh: &TetMesh<V, T>) -> VertexId {
        unsafe { mesh.tets.get_unchecked(self.tet).vertices[self.edge as usize % 4] }
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

    /// Iterates over the walkers with the same current edge as this one
    /// (including this one) in ring order.
    pub fn edge_ring<'a, V, T>(self, mesh: &'a TetMesh<V, T>) -> WalkerEdgeRing<'a, V, T> {
        WalkerEdgeRing(CircularListIter::new(
            mesh,
            self,
            |_, _| (),
            |mesh, w| w.to_twin_edge().to_adj(mesh),
        ))
    }

    /// Gets the current tet id.
    pub fn id(self) -> TetId {
        self.tet
    }

    /// Gets the current undirected edge index.
    fn undir_edge_index(self) -> usize {
        unsafe { *[0, 0, 1, 1, 2, 3, 3, 2, 4, 5, 4, 5].get_unchecked(self.edge as usize) }
    }

    /// Gets the current undirected edge id.
    fn undir_edge<V, T>(self, mesh: &TetMesh<V, T>) -> UndirEdgeId {
        unsafe { mesh.tets.get_unchecked(self.id()).edges[self.undir_edge_index()] }
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
            unsafe {
                self.edge = *[10, 11, 8, 9, 1, 0, 3, 2, 7, 6, 5, 4].get_unchecked(self.edge as usize);
            }
            self
        }
    }

    crate::alias! {
        to_pve,
        /// Advance to the previous edge from the current vertex. This is a clockwise rotation.
        pub fn to_prev_vertex_edge((mut) self: Self) -> Self {
            unsafe {
                self.edge = *[5, 4, 7, 6, 11, 10, 9, 8, 2, 3, 0, 1].get_unchecked(self.edge as usize);
            }
            self
        }
    }

    /// Flip the current edge, keeping the current tet the same.
    pub fn to_twin_edge(mut self) -> Self {
        unsafe {
            self.edge = *[1, 0, 3, 2, 7, 6, 5, 4, 10, 11, 8, 9].get_unchecked(self.edge as usize);
        }
        self
    }

    /// Sets the current edge to the opposite edge of the current tet.
    pub fn to_opp_edge(mut self) -> Self {
        unsafe {
            self.edge = *[2, 3, 0, 1, 5, 4, 7, 6, 11, 10, 9, 8].get_unchecked(self.edge as usize);
        }
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
        unsafe {
            // Don't question the indexing order here, it just works™
            self.edge = Self::PERMUTATION_FWD.get_unchecked(self.edge as usize)[perm as usize];
        }
        self
    }

    crate::alias! {
        to_adj,
        /// Flip the current triangle, moving to the adjacent tet on that triangle.
        /// This flips the current edge.
        pub fn to_twin_tri<V, T>((mut) self: Self, mesh: &TetMesh<V, T>) -> Self {
            let div = self.edge / 4;
            unsafe {
                self = mesh.tets.get_unchecked(self.id()).opp_tets[self.edge as usize % 4];
            }
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
            unsafe {
                mesh.tets.get_unchecked(self.id()).opp_tets[self.edge as usize % 4]
            }
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

    fn to_edge(mut self, edge: u32) -> Self {
        self.edge = edge;
        self
    }

    /// Walk to a specific edge assuming the walker is on the triangle containing that edge.
    fn to_edge_assuming_tri<V, T>(self, mesh: &TetMesh<V, T>, edge: [VertexId; 2]) -> Self {
        let idx = self.tet(mesh).iter().position(|v| *v == edge[0]).unwrap();
        match idx {
            0 => self,
            1 => self.to_nfe(),
            2 => self.to_pfe(),
            _ => unreachable!(),
        }
    }

    /// Assumes this walker is on a sliver tet and gets the kind of sliver it's on.
    fn sliver_kind<V, T>(self, mesh: &mut TetMesh<V, T>) -> SliverKind {
        let p0 = mesh.vertices[self.first(mesh)].position();
        let p1 = mesh.vertices[self.second(mesh)].position();
        let p2 = mesh.vertices[self.third(mesh)].position();
        let p3 = mesh.vertices[self.fourth(mesh)].position();

        let t0 = (p2 - p3).cross(&(p1 - p3)).norm();
        let t1 = (p3 - p2).cross(&(p0 - p2)).norm();
        let t2 = (p0 - p1).cross(&(p3 - p1)).norm();
        let t3 = (p1 - p0).cross(&(p2 - p0)).norm();
        let e03 = (p0 - p3).cross(&(p1 - p2)).norm();
        let e13 = (p1 - p3).cross(&(p2 - p0)).norm();
        let e23 = (p2 - p3).cross(&(p0 - p1)).norm();
        let max = t0.max(t1).max(t2).max(t3).max(e03).max(e13).max(e23);

        match max {
            _ if max == t0 => SliverKind::TriVertex(self.to_perm(Permutation::_3210)),
            _ if max == t1 => SliverKind::TriVertex(self.to_opp_edge()),
            _ if max == t2 => SliverKind::TriVertex(self.to_twin_edge()),
            _ if max == t3 => SliverKind::TriVertex(self),
            _ if max == e03 => SliverKind::EdgeEdge(self.to_pve()),
            _ if max == e13 => SliverKind::EdgeEdge(self.to_nve()),
            _ if max == e23 => SliverKind::EdgeEdge(self),
            _ => unreachable!()
        }
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

        if mesh.track_tri_map {
            for i in 0..4 {
                mesh.tri_map
                    .remove(&TriVertices::new(self.to_edge(i).tri(mesh)));
            }
        }

        // The vertex didn't have a tet before.
        mesh.vertices[vertex].tet = self;

        // Just in case of dangling vertex-tet reference
        let vs = mesh.tets[id0].vertices();
        mesh.vertices
            .get_mut(vs[0])
            .map(|v| v.tet = TetWalker::new(id2, 3));
        mesh.vertices
            .get_mut(vs[1])
            .map(|v| v.tet = TetWalker::new(id3, 2));
        mesh.vertices
            .get_mut(vs[2])
            .map(|v| v.tet = TetWalker::new(id0, 1));
        mesh.vertices
            .get_mut(vs[3])
            .map(|v| v.tet = TetWalker::new(id1, 0));

        // Add new edges
        let es = mesh
            .edges
            .extend_values((0..4).map(|i| UndirEdge::new([vs[i], vertex])))
            .collect::<Vec<_>>();

        // Update tets
        // Note that adjacent tets need to change their opposite tet ids,
        // but not their opposite tet edge indexes because they refer to the same face.
        mesh.tets[id0].vertices[0] = vertex;
        mesh.tets[id0].opp_tets[1] = TetWalker::new(id1, 0);
        mesh.tets[id0].opp_tets[2] = TetWalker::new(id2, 4);
        mesh.tets[id0].opp_tets[3] = TetWalker::new(id3, 8);
        mesh.tets[id0].edges[1] = es[1];
        mesh.tets[id0].edges[5] = es[2];
        mesh.tets[id0].edges[3] = es[3];

        mesh.tets[id1].vertices[1] = vertex;
        mesh.tets[id1].opp_tets[2] = TetWalker::new(id2, 9);
        mesh.tets[id1].opp_tets[3] = TetWalker::new(id3, 5);
        mesh.tets[id1].opp_tets[0] = TetWalker::new(id0, 1);
        mesh.tets[id1].edges[2] = es[2];
        mesh.tets[id1].edges[4] = es[3];
        mesh.tets[id1].edges[1] = es[0];

        mesh.tets[id2].vertices[2] = vertex;
        mesh.tets[id2].opp_tets[3] = TetWalker::new(id3, 2);
        mesh.tets[id2].opp_tets[0] = TetWalker::new(id0, 6);
        mesh.tets[id2].opp_tets[1] = TetWalker::new(id1, 10);
        mesh.tets[id2].edges[0] = es[3];
        mesh.tets[id2].edges[5] = es[0];
        mesh.tets[id2].edges[2] = es[1];

        mesh.tets[id3].vertices[3] = vertex;
        mesh.tets[id3].opp_tets[0] = TetWalker::new(id0, 11);
        mesh.tets[id3].opp_tets[1] = TetWalker::new(id1, 7);
        mesh.tets[id3].opp_tets[2] = TetWalker::new(id2, 3);
        mesh.tets[id3].edges[3] = es[0];
        mesh.tets[id3].edges[4] = es[1];
        mesh.tets[id3].edges[0] = es[2];

        // Update adjacent tets' opposite indexes
        let ids = [id0, id1, id2, id3];
        for i in 1..4 {
            let walker = mesh[ids[i]].opp_tets[i];
            mesh.tets[walker.id()].opp_tets[walker.edge as usize % 4].tet = ids[i];
        }

        if mesh.track_tri_map {
            for id in &ids {
                for i in 0..4 {
                    let walker = TetWalker::new(*id, i);
                    mesh.tri_map
                        .insert(TriVertices::new(walker.tri(mesh)), walker);
                }
            }
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
    /// Returns walkers of the 3 tets created. Each walker's opposite edge is the new edge
    /// and the first walker's current edge is the same as this walker's current edge.
    pub fn flip23_unchecked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> [TetWalker; 3] {
        // Walker on other tri
        let other = self.to_adj(&mesh);

        let id0 = self.id();
        let id1 = other.id();
        let id2 = mesh.tets.insert(mesh.tets[id0].clone());
        let ids = [id0, id1, id2];

        if mesh.track_tri_map {
            for walker in &[self, other] {
                for i in 0..4 {
                    mesh.tri_map
                        .remove(&TriVertices::new(walker.to_edge(i).tri(mesh)));
                }
            }
        }

        // Fix tet walkers on adjacent tets
        let mut adj_t = [TetWalker::new(TetId::invalid(), 0); 3];
        for (i, perm) in [Permutation::_1032, Permutation::_0231, Permutation::_2130]
            .iter()
            .enumerate()
        {
            let walk_t = self.to_perm(*perm);
            adj_t[i] = walk_t.to_adj(&mesh);
            let walker = &mut mesh.tets[adj_t[i].id()].opp_tets[adj_t[i].edge as usize % 4];
            *walker = TetWalker::new(
                ids[i],
                Self::PERMUTATION_INV[walk_t.edge as usize][walker.edge as usize],
            );
        }

        // Bottom tets; permutations are a little different
        let mut adj_b = [TetWalker::new(TetId::invalid(), 0); 3];
        for (i, perm) in [Permutation::_1032, Permutation::_2130, Permutation::_0231]
            .iter()
            .enumerate()
        {
            adj_b[i] = other.to_perm(*perm).to_adj(&mesh);
            let walker = &mut mesh.tets[adj_b[i].id()].opp_tets[adj_b[i].edge as usize % 4];
            let edge = match i {
                0 => other.edge,
                1 => other.to_nfe().edge,
                2 => other.to_pfe().edge,
                _ => unsafe { unreachable_unchecked() }, // i is in 0..3
            };
            *walker = TetWalker::new(
                ids[i],
                Self::PERMUTATION_INV[edge as usize][walker.edge as usize],
            );
        }

        // Introduce the 3 new tets
        let ([v1, v0, v2, v3], v4) = (self.tet(&mesh), other.fourth(&mesh));

        // Just in case of dangling vertex-tet reference
        mesh.vertices
            .get_mut(v0)
            .map(|v| v.tet = TetWalker::new(id0, 3));
        mesh.vertices
            .get_mut(v1)
            .map(|v| v.tet = TetWalker::new(id1, 3));
        mesh.vertices
            .get_mut(v2)
            .map(|v| v.tet = TetWalker::new(id2, 3));
        mesh.vertices
            .get_mut(v3)
            .map(|v| v.tet = TetWalker::new(id0, 1));
        mesh.vertices
            .get_mut(v4)
            .map(|v| v.tet = TetWalker::new(id0, 0));

        let edge = mesh.edges.insert(UndirEdge::new([v3, v4]));

        for (i, [va, vb]) in [[v0, v1], [v1, v2], [v2, v0]].iter().enumerate() {
            let tet = &mut mesh.tets[ids[i]];
            tet.vertices = [*va, *vb, v3, v4];
            tet.opp_tets[0] = TetWalker::new(ids[(i + 1) % 3], 1);
            tet.opp_tets[1] = TetWalker::new(ids[(i + 2) % 3], 0);
            tet.opp_tets[2] = TetWalker::new(adj_b[i].id(), adj_b[i].edge);
            tet.opp_tets[3] = TetWalker::new(adj_t[i].id(), adj_t[i].edge);
            tet.edges[0] = edge;

            // Fix undirected edge ids
            mesh.tets[ids[i]].edges[1] = adj_t[i].undir_edge(mesh);
            mesh.tets[ids[i]].edges[2] = adj_t[i].to_pfe().undir_edge(mesh);
            mesh.tets[ids[i]].edges[3] = adj_b[i].to_pfe().undir_edge(mesh);
            mesh.tets[ids[i]].edges[4] = adj_b[i].to_nfe().undir_edge(mesh);
            mesh.tets[ids[i]].edges[5] = adj_t[i].to_nfe().undir_edge(mesh);
        }

        if mesh.track_tri_map {
            for id in &ids {
                for i in 0..4 {
                    let walker = TetWalker::new(*id, i);
                    mesh.tri_map
                        .insert(TriVertices::new(walker.tri(mesh)), walker);
                }
            }
        }

        [
            TetWalker::new(id0, 2),
            TetWalker::new(id1, 2),
            TetWalker::new(id2, 2),
        ]
    }

    /// Performs a 2-to-3 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets, but does not flip away locked tris.
    /// Returns walkers of the 3 tets created. Each walker's opposite edge is the new edge
    /// and the first walker's current edge is the same as this walker's current edge.
    pub fn flip23_unlocked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> Option<[TetWalker; 3]> {
        let edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::LOCKED_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::LOCKED_0.bits << edge % 4)
        };

        if mesh[self.id()].flags.get().intersects(edge_flag(self.edge)) {
            None
        } else {
            Some(self.flip23_unchecked(mesh))
        }
    }

    /// Performs a 2-to-3 flip on the mesh, but only if doing so improves the quality.
    /// Returns walkers of the 3 tets created. Each walker's opposite edge is the new edge
    /// and the first walker's current edge is the same as this walker's current edge.
    pub fn flip23_quality<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> Option<[TetWalker; 3]> {
        let other = self.to_adj_ae(mesh);
        let ([v1, v0, v2, v3], v4) = (self.tet(&mesh), other.fourth(&mesh));

        let q0 = mesh.quality(v1, v0, v2, v3);
        let q1 = mesh.quality(v0, v1, v2, v4);
        let r0 = mesh.quality(v0, v1, v3, v4);
        let r1 = mesh.quality(v1, v2, v3, v4);
        let r2 = mesh.quality(v2, v0, v3, v4);

        if r0.min(r1).min(r2) > q0.min(q1) {
            self.flip23_unlocked(mesh)
        } else {
            None
        }
    }

    /// Performs a 3-to-2 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Returns walkers for the 2 tets created by the flip. Each walker's opposite triangle is the new triangle.
    pub fn flip32_unchecked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> [TetWalker; 2] {
        let other1 = self.to_twin_edge().to_adj(&mesh);
        let other2 = other1.to_twin_edge().to_adj(&mesh);
        let old_walkers = [self, other1, other2];

        let id0 = self.id();
        let id1 = other1.id();

        if mesh.track_tri_map {
            for walker in &[self, other1, other2] {
                for i in 0..4 {
                    mesh.tri_map
                        .remove(&TriVertices::new(walker.to_edge(i).tri(mesh)));
                }
            }
        }

        // Fix tet walkers on adjacent tets
        let mut adj_t = [TetWalker::new(TetId::invalid(), 0); 3];
        for (i, perm) in [Permutation::_3210, Permutation::_2130, Permutation::_1320]
            .iter()
            .enumerate()
        {
            adj_t[i] = old_walkers[i].to_opp_edge().to_adj(&mesh);
            let walker = &mut mesh.tets[adj_t[i].id()].opp_tets[adj_t[i].edge as usize % 4];
            *walker = TetWalker::new(
                id0,
                Self::PERMUTATION_INV[old_walkers[i].to_perm(*perm).edge as usize]
                    [walker.edge as usize],
            );
        }

        // Bottom
        let mut adj_b = [TetWalker::new(TetId::invalid(), 0); 3];
        for (i, perm) in [Permutation::_2301, Permutation::_0231, Permutation::_3021]
            .iter()
            .enumerate()
        {
            adj_b[i] = old_walkers[i].to_perm(Permutation::_3210).to_adj(&mesh);
            let walker = &mut mesh.tets[adj_b[i].id()].opp_tets[adj_b[i].edge as usize % 4];
            *walker = TetWalker::new(
                id1,
                Self::PERMUTATION_INV[old_walkers[i].to_perm(*perm).edge as usize]
                    [walker.edge as usize],
            );
        }

        // Introduce the 2 new tets
        let ([v3, v4, v1, v0], v2) = (self.tet(&mesh), other1.fourth(&mesh));

        // Just in case of dangling vertex-tet reference
        mesh.vertices
            .get_mut(v0)
            .map(|v| v.tet = TetWalker::new(id0, 3));
        mesh.vertices
            .get_mut(v1)
            .map(|v| v.tet = TetWalker::new(id0, 7));
        mesh.vertices
            .get_mut(v2)
            .map(|v| v.tet = TetWalker::new(id0, 11));
        mesh.vertices
            .get_mut(v3)
            .map(|v| v.tet = TetWalker::new(id0, 0));
        mesh.vertices
            .get_mut(v4)
            .map(|v| v.tet = TetWalker::new(id1, 0));

        let tet = &mut mesh.tets[id0];
        tet.vertices = [v0, v1, v2, v3];
        tet.opp_tets[0] = TetWalker::new(adj_t[2].id(), adj_t[2].to_nfe().edge);
        tet.opp_tets[1] = TetWalker::new(adj_t[1].id(), adj_t[1].to_pfe().edge);
        tet.opp_tets[2] = TetWalker::new(adj_t[0].id(), adj_t[0].edge);
        tet.opp_tets[3] = TetWalker::new(id1, 3);
        mesh.tets[id0].edges[0] = adj_t[1].to_pfe().undir_edge(mesh);
        mesh.tets[id0].edges[1] = adj_t[0].undir_edge(mesh);
        mesh.tets[id0].edges[2] = adj_t[2].undir_edge(mesh);
        mesh.tets[id0].edges[3] = adj_t[0].to_pfe().undir_edge(mesh);
        mesh.tets[id0].edges[4] = adj_t[2].to_pfe().undir_edge(mesh);
        mesh.tets[id0].edges[5] = adj_t[1].undir_edge(mesh);

        let tet = &mut mesh.tets[id1];
        tet.vertices = [v1, v0, v2, v4];
        tet.opp_tets[0] = TetWalker::new(adj_b[1].id(), adj_b[1].to_nfe().edge);
        tet.opp_tets[1] = TetWalker::new(adj_b[2].id(), adj_b[2].to_pfe().edge);
        tet.opp_tets[2] = TetWalker::new(adj_b[0].id(), adj_b[0].edge);
        tet.opp_tets[3] = TetWalker::new(id0, 3);
        mesh.tets[id1].edges[0] = adj_b[1].to_nfe().undir_edge(mesh);
        mesh.tets[id1].edges[1] = mesh.tets[id0].edges[1];
        mesh.tets[id1].edges[2] = mesh.tets[id0].edges[5];
        mesh.tets[id1].edges[3] = adj_b[2].to_nfe().undir_edge(mesh);
        mesh.tets[id1].edges[4] = adj_b[0].to_nfe().undir_edge(mesh);
        mesh.tets[id1].edges[5] = mesh.tets[id0].edges[2];

        mesh.edges.remove(other2.undir_edge(mesh));
        mesh.tets.remove(other2.id());

        if mesh.track_tri_map {
            for id in &[id0, id1] {
                for i in 0..4 {
                    let walker = TetWalker::new(*id, i);
                    mesh.tri_map
                        .insert(TriVertices::new(walker.tri(mesh)), walker);
                }
            }
        }

        [TetWalker::new(id0, 10), TetWalker::new(id1, 10)]
    }

    /// Performs a 3-to-2 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets, but does not flip locked edges or tris.
    /// This also checks that there are only 3 tets around the edge to flip away.
    /// Returns walkers for the 2 tets created by the flip. Each walker's opposite triangle is the new triangle.
    pub fn flip32_unlocked<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> Option<[TetWalker; 2]> {
        let edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::LOCKED_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::LOCKED_0.bits << edge % 4)
        };

        let other1 = self.to_twin_edge().to_adj(&mesh);
        let other2 = other1.to_twin_edge().to_adj(&mesh);
        if mesh[self.id()].flags.get().intersects(edge_flag(self.edge)) ||
            mesh[other1.id()].flags.get().intersects(edge_flag(other1.edge)) ||
            mesh[other2.id()].flags.get().intersects(edge_flag(other2.edge)) ||
            mesh.edges[self.undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED) ||
            self.to_twin_edge().to_adj(mesh).to_twin_edge().to_adj(mesh).to_twin_edge().to_adj(mesh) != self
        {
            None
        } else {
            Some(self.flip32_unchecked(mesh))
        }
    }

    /// Performs a 3-to-2 flip on the mesh, but only if it improves the quality.
    /// Returns walkers for the 2 tets created by the flip. Each walker's opposite triangle is the new triangle.
    pub fn flip32_quality<V: Clone, T: Clone>(self, mesh: &mut TetMesh<V, T>) -> Option<[TetWalker; 2]> {
        let other1 = self.to_twin_edge().to_adj(&mesh);
        let ([v3, v4, v1, v0], v2) = (self.tet(&mesh), other1.fourth(&mesh));

        let q0 = mesh.quality(v3, v4, v1, v0);
        let q1 = mesh.quality(v3, v4, v0, v2);
        let q2 = mesh.quality(v3, v4, v2, v1);
        let r0 = mesh.quality(v0, v1, v2, v3);
        let r1 = mesh.quality(v1, v0, v2, v4);

        if r0.min(r1) > q0.min(q1).min(q2) {
            self.flip32_unlocked(mesh)
        } else {
            None
        }
    }

    /// Attempts to remove this walker's edge to make progress towards removing something from `goals`.
    /// Returns the index in `goals` to the removed edge and a walker to the new tri if it succeeds or something from `goals` is removed.
    /// The vertices in `goals` must be sorted.
    /// This function requires that the mesh is tracking the triangle map.
    ///
    /// - `bad_edges` is a set of edges to not add tris to.
    ///
    /// WARNING: This invalidates all tet walkers except those tracked by `mesh.tri_map`.
    fn progress_by_removing_edge<
        V: Clone,
        T: Clone,
        EA: Fn(&TetMesh<V, T>, [VertexId; 2]) -> bool + Clone,
        FA: Fn(&TetMesh<V, T>, [VertexId; 3]) -> bool + Clone,
    >(
        self,
        mesh: &mut TetMesh<V, T>,
        edge_addable: EA,
        tri_addable: FA,
        depth: usize,
        allow_slivers: bool,
        goals: &mut Vec<EdgeOrTri>,
        bad_edges: &mut Vec<[VertexId; 2]>,
    ) -> Option<usize> {
        if depth == 0
            || mesh.edges[self.undir_edge(mesh)]
                .flags
                .get()
                .intersects(UndirEdgeFlags::LOCKED)
        {
            return None;
        }

        // Be sure to pop this goal before every return
        goals.push(EdgeOrTri::Edge(sorted_2(self.edge(mesh))));

        // Don't use the walkers directly since they get invalidated after each flip.
        let edge = self.edge(mesh);
        let mut tri_opps = self
                .edge_ring(mesh)
                .map(|w| w.third(mesh))
                .collect::<Vec<_>>();
        let mut next_tri_opps = vec![];

        // Can't do a 3-to-2 flip without starting from 3
        'reduce: while {
            tri_opps.append(&mut next_tri_opps);
            tri_opps.len() > 3
        } {
            let initial = tri_opps.len();

            // See if we can reduce the number of triangles around the edge by ≥1
            while let Some(v2) = tri_opps.pop() {
                if let Some(walker) = mesh.tri_map.get(&TriVertices::new([edge[0], edge[1], v2])) {
                    let walker = walker.to_edge_assuming_tri(mesh, edge);

                    // Avoid flipping a triangle to the edge I'm trying to remove triangles from!
                    bad_edges.push(sorted_2(edge));
                    let index = walker.progress_by_removing_tri(
                        mesh,
                        edge_addable.clone(),
                        tri_addable.clone(),
                        depth - 1,
                        allow_slivers,
                        goals,
                        bad_edges,
                    );
                    bad_edges.pop();

                    // Failure to remove triangle means try the next one
                    if index.map(|i| i < goals.len()).unwrap_or(false) {
                        goals.pop();
                        return index;
                    } else if index.is_none() {
                        next_tri_opps.push(v2);
                    }
                } 

                if tri_opps.len() + next_tri_opps.len() <= 3 {
                    // Make sure to obtain the new tri_opps
                    continue 'reduce;
                }
            }

            if initial == next_tri_opps.len() {
                goals.pop();
                return None;
            }
        }
        assert_eq!(tri_opps.len(), 3);

        let edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::LOCKED_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::LOCKED_0.bits << edge % 4)
        };

        // Edge hasn't been removed yet, and there are 3 triangles around the edge, so this should work.
        let walk_0 = mesh.tri_map[&TriVertices::new([edge[0], edge[1], *tri_opps.last().unwrap()])]
            .to_edge_assuming_tri(mesh, edge);
        let walk_1 = walk_0.to_twin_edge().to_adj(mesh);
        let walk_2 = walk_1.to_twin_edge().to_adj(mesh);

        let [v3, v4] = edge;
        let ([v0, v1], v2) = (walk_0.opp_edge(mesh), walk_1.fourth(mesh));

        // Concavity => try to remove a triangle.
        // If it can be removed, so can the edge.
        let test = if allow_slivers {
            !mesh.orient_3d(v0, v1, v3, v2) || !mesh.orient_3d(v0, v1, v2, v4)
        } else {
            !mesh.non_sliver(v0, v1, v3, v2) || !mesh.non_sliver(v0, v1, v2, v4)
        };
        if test {
            for v in &[v0, v1, v2] {
                // Walker should exist if previous flip didn't work.
                let walker = mesh.tri_map[&TriVertices::new([edge[0], edge[1], *v])]
                    .to_edge_assuming_tri(mesh, edge);

                // Avoid flipping a triangle to the edge I'm trying to remove triangles from!
                bad_edges.push(sorted_2(edge));
                let index = walker.progress_by_removing_tri(
                    mesh,
                    edge_addable.clone(),
                    tri_addable.clone(),
                    depth - 1,
                    allow_slivers,
                    goals,
                    bad_edges,
                );
                bad_edges.pop();

                // Failure to remove triangle means try the next one
                if index.is_some() {
                    goals.pop();
                    return index;
                }
            }
            goals.pop();
            return None;
        }

        // Check removability
        if !bad_edges.contains(&sorted_2([v0, v1]))
            && !bad_edges.contains(&sorted_2([v1, v2]))
            && !bad_edges.contains(&sorted_2([v2, v0]))
            && tri_addable(mesh, [v0, v1, v2])
            && !mesh.tets[walk_0.id()]
                .flags
                .get()
                .intersects(edge_flag(walk_0.edge))
            && !mesh.tets[walk_1.id()]
                .flags
                .get()
                .intersects(edge_flag(walk_1.edge))
            && !mesh.tets[walk_2.id()]
                .flags
                .get()
                .intersects(edge_flag(walk_2.edge))
        {
            walk_0.flip32_unchecked(mesh);

            // Get the *first* index of the goal, in case an earlier function call wanted to remove the same edge.
            // If the goal is not here, that's a bug.
            let goal_index = goals
                .iter()
                .position(|g| {
                    [
                        *goals.last().unwrap(),
                        EdgeOrTri::Tri(sorted_3([v0, v3, v4])),
                        EdgeOrTri::Tri(sorted_3([v1, v3, v4])),
                        EdgeOrTri::Tri(sorted_3([v2, v3, v4])),
                    ]
                    .contains(g)
                })
                .unwrap();
            goals.pop();
            Some(goal_index)
        } else {
            goals.pop();
            None
        }
    }

    /// Attempts to remove this walker's triangle to make progress towards removing something from `goals`.
    /// Returns the index in `goals` to the removed triangle and a walker to the same edge if it succeeds or something from `goals` is removed.
    /// The vertices in `goals` must be sorted.
    /// This function requires that the mesh is tracking the triangle map.
    ///
    /// - `bad_edges` is a set of edges to not add tris to.
    ///
    /// WARNING: This invalidates all tet walkers except those tracked by `mesh.tri_map`.
    fn progress_by_removing_tri<
        V: Clone,
        T: Clone,
        EA: Fn(&TetMesh<V, T>, [VertexId; 2]) -> bool + Clone,
        FA: Fn(&TetMesh<V, T>, [VertexId; 3]) -> bool + Clone,
    >(
        self,
        mesh: &mut TetMesh<V, T>,
        edge_addable: EA,
        tri_addable: FA,
        depth: usize,
        allow_slivers: bool,
        goals: &mut Vec<EdgeOrTri>,
        bad_edges: &mut Vec<[VertexId; 2]>,
    ) -> Option<usize> {
        if depth == 0 {
            return None;
        }

        let edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::LOCKED_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::LOCKED_0.bits << edge % 4)
        };
        if mesh.tets[self.id()]
            .flags
            .get()
            .intersects(edge_flag(self.edge))
        {
            return None;
        }

        // Be sure to pop this goal before every return
        goals.push(EdgeOrTri::Tri(sorted_3(self.tri(mesh))));

        // Make sure there's no edge blocking the way.
        let adj = self.to_adj_ae(mesh);
        let ([v0, v1, v2, v3], v4) = (self.tet(mesh), adj.fourth(mesh));

        for (va, vb, walker) in &[
            (v0, v1, self),
            (v1, v2, self.to_nfe()),
            (v2, v0, self.to_pfe()),
        ] {
            let test = if allow_slivers {
                !mesh.orient_3d(*va, *vb, v4, v3)
            } else {
                !mesh.non_sliver(*va, *vb, v4, v3)
            };
            if test {
                let index = walker.progress_by_removing_edge(
                    mesh,
                    edge_addable,
                    tri_addable,
                    depth - 1,
                    allow_slivers,
                    goals,
                    bad_edges,
                );

                goals.pop();
                return index;
            }
        }

        if !bad_edges.contains(&sorted_2([v0, v3]))
            && !bad_edges.contains(&sorted_2([v0, v4]))
            && !bad_edges.contains(&sorted_2([v1, v3]))
            && !bad_edges.contains(&sorted_2([v1, v4]))
            && !bad_edges.contains(&sorted_2([v2, v3]))
            && !bad_edges.contains(&sorted_2([v2, v4]))
            && tri_addable(mesh, [v0, v3, v4])
            && tri_addable(mesh, [v1, v3, v4])
            && tri_addable(mesh, [v2, v3, v4])
            && edge_addable(mesh, [v3, v4])
        {
            // No flip happened yet, so using `self` is fine.
            self.flip23_unchecked(mesh);

            // Get the *first* index of the goal, in case an earlier function call wanted to remove the same tri.
            // If the goal is not here, that's a bug.
            let goal_index = goals
                .iter()
                .position(|g| g == goals.last().unwrap())
                .unwrap();
            goals.pop();
            Some(goal_index)
        } else {
            goals.pop();
            None
        }
    }

    /// Returns the boundary of the region of tets that will be deleted by removing some vertex, starting at the current tet.
    /// Also returns all tets in the region.
    /// Each item in the boundary is a tet walker whose current triangle is a triangle of the boundary.
    pub fn point_cavity<V, T>(
        self,
        mesh: &TetMesh<V, T>,
        convex: bool,
        vertex: VertexId,
    ) -> (Vec<TetWalker>, Vec<TetId>) {
        // Initialize intermediate boundary
        let mut bound_imm = vec![
            TetWalker::new(self.id(), 0),
            TetWalker::new(self.id(), 1),
            TetWalker::new(self.id(), 2),
            TetWalker::new(self.id(), 3),
        ].into_iter().collect::<VecDeque<_>>();
        mesh.tets[self.id()].set_flags(
            TetFlags::BOUND_IMM_0
                | TetFlags::BOUND_IMM_1
                | TetFlags::BOUND_IMM_2
                | TetFlags::BOUND_IMM_3,
        );

        let mut boundary = vec![];
        let mut enclosed = vec![self.id()];

        let imm_edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::BOUND_IMM_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::BOUND_IMM_0.bits << edge % 4)
        };
        let limit_edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::BOUND_LIMIT_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::BOUND_LIMIT_0.bits << edge % 4)
        };
        let locked_edge_flag = |edge| unsafe {
            // Safety: tri.edge % 4 is in 0..4, and TetFlags::LOCKED_0.bits << i for i in 0..4 is a valid flag.
            TetFlags::from_bits_unchecked(TetFlags::LOCKED_0.bits << edge % 4)
        };

        // Extend intermediate boundary
        while let Some(tri) = bound_imm.pop_front() {
            // Clear flag
            let tet = &mesh.tets[tri.id()];
            if !tet.flags.get().intersects(imm_edge_flag(tri.edge)) {
                // Triangle was removed by twin triangle.
                continue;
            }
            tet.clear_flags(imm_edge_flag(tri.edge));

            let adj = tri.to_adj_ae(mesh);
            let vs = adj.tet(mesh);

            // Check that no edge that will be removed is locked
            let removal_safe = if convex {
                true
            } else {
                let walk0 = adj.to_perm(Permutation::_3210).to_adj_ae(mesh);
                let flags0 = &mesh.tets[walk0.id()].flags;
                let walk1 = adj.to_opp_edge().to_adj_ae(mesh);
                let flags1 = &mesh.tets[walk1.id()].flags;
                let walk2 = adj.to_twin_edge().to_adj_ae(mesh);
                let flags2 = &mesh.tets[walk2.id()].flags;

                let rm0 = flags0.get().intersects(imm_edge_flag(walk0.edge));
                let rm1 = flags1.get().intersects(imm_edge_flag(walk1.edge));
                let rm2 = flags2.get().intersects(imm_edge_flag(walk2.edge));
                
                // Don't extend past the limit
                !flags0.get().intersects(limit_edge_flag(walk0.edge)) &&
                !flags1.get().intersects(limit_edge_flag(walk1.edge)) &&
                !flags2.get().intersects(limit_edge_flag(walk2.edge)) &&
                // Make sure edges are removable
                (!rm0 || !mesh.edges[adj.to_perm(Permutation::_1203).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                (!rm1 || !mesh.edges[adj.to_perm(Permutation::_2013).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                (!rm2 || !mesh.edges[adj.to_perm(Permutation::_0123).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                (!rm0 || !rm1 || !mesh.edges[adj.to_perm(Permutation::_2301).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                (!rm1 || !rm2 || !mesh.edges[adj.to_perm(Permutation::_0312).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                (!rm2 || !rm0 || !mesh.edges[adj.to_perm(Permutation::_1320).undir_edge(mesh)].flags.get().intersects(UndirEdgeFlags::LOCKED)) &&
                // Make sure tri is removable
                !mesh[adj.id()].flags.get().intersects(locked_edge_flag(adj.edge))
            };

            if removal_safe && (mesh.in_sphere_sim(vs[0], vs[1], vs[2], vs[3], vertex) ||
                (!convex && mesh.tet_close_to_vertex(vs[0], vs[1], vs[2], vs[3], vertex))) {
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
                    let twin = walker.to_adj_ae(mesh);
                    if mesh.tets[twin.id()]
                        .flags
                        .get()
                        .intersects(imm_edge_flag(twin.edge))
                    {
                        mesh.tets[twin.id()].clear_flags(imm_edge_flag(twin.edge));
                    } else {
                        mesh.tets[walker.id()].set_flags(imm_edge_flag(walker.edge));
                        bound_imm.push_back(*walker);
                    }
                }
            } else {
                // Reached edge of boundary
                boundary.push(tri);
                if !convex {
                    // Mark so the expansion doesn't try to gobble this up.
                    mesh.tets[tri.id()].set_flags(limit_edge_flag(tri.edge));
                }
            }
        }

        if !convex {
            // Region might not be star-shaped, so fix that.
            for walker in boundary.drain(..) {
                bound_imm.push_back(walker);
                mesh.tets[walker.id()].set_flags(imm_edge_flag(walker.edge));
                mesh.tets[walker.id()].clear_flags(limit_edge_flag(walker.edge));
            }

            while let Some(tri) = bound_imm.pop_front() {
                if !mesh.tets[tri.id()].flags.get().intersects(imm_edge_flag(tri.edge)) {
                    // Triangle was removed by twin triangle.
                    continue;
                }

                let vs = tri.tri(mesh);
                if mesh.orient_3d(vs[0], vs[1], vs[2], vertex) {
                    boundary.push(tri);
                } else {
                    // Remove from cavity
                    let tet = &mesh.tets[tri.id()];
                    tet.set_flags(TetFlags::NOT_ENCLOSURE);
                    tet.clear_flags(imm_edge_flag(tri.edge));
                    
                    // Local retraction
                    for walker in [
                        tri.to_opp_edge(),
                        tri.to_twin_edge(),
                        tri.to_perm(Permutation::_3210),
                    ]
                    .iter()
                    {
                        if mesh.tets[walker.id()]
                            .flags
                            .get()
                            .intersects(imm_edge_flag(walker.edge))
                        {
                            mesh.tets[walker.id()].clear_flags(imm_edge_flag(walker.edge));
                        } else {
                            let twin = walker.to_adj_ae(mesh);
                            mesh.tets[twin.id()].set_flags(imm_edge_flag(twin.edge));
                            bound_imm.push_back(twin);
                        }
                    }
                }
            }

            // Update boundary and enclosure.
            boundary = boundary.into_iter().filter(|walker| {
                let flag = mesh.tets[walker.id()].flags.get().intersects(imm_edge_flag(walker.edge));
                mesh.tets[walker.id()].clear_flags(imm_edge_flag(walker.edge));
                flag
            }).collect();

            enclosed = enclosed.into_iter().filter(|id| {
                let flag = mesh.tets[*id].flags.get().intersects(TetFlags::NOT_ENCLOSURE);
                mesh.tets[*id].clear_flags(TetFlags::NOT_ENCLOSURE);
                !flag
            }).collect();

            assert!(!enclosed.is_empty(), "Adding vertex {} removed 0 tets.", vertex);
        }

        debug_assert!(
            enclosed.len() == enclosed.iter().collect::<FnvHashSet<_>>().len(),
            "Enclosure contains repeat elements unexpectedly."
        );
        (boundary, enclosed)
    }
}

/// Iterates over tet walkers with the same current edge in ring order.
#[derive(Clone, Debug)]
pub struct WalkerEdgeRing<'a, V, T>(
    CircularListIter<
        &'a TetMesh<V, T>,
        TetWalker,
        fn(&'a TetMesh<V, T>, TetWalker) -> (),
        fn(&'a TetMesh<V, T>, TetWalker) -> TetWalker,
    >,
);

impl<'a, V, T> Iterator for WalkerEdgeRing<'a, V, T> {
    type Item = TetWalker;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|(walker, _)| walker)
    }
}

/// A manifold tetrahedralization.
#[derive(Clone, Debug)]
pub struct TetMesh<V, T> {
    vertices: IdMap<VertexId, Vertex<V>>,
    ghost_flags: Cell<VertexFlags>,
    edges: IdMap<UndirEdgeId, UndirEdge>,
    tets: IdMap<TetId, Tet<T>>,
    /// Map from even-sorted triangles to tet walkers for those triangles.
    /// This is slow to track, so don't track it all the time.
    tri_map: FnvHashMap<TriVertices, TetWalker>,
    track_tri_map: bool,
    /// Relevant after the initial Delaunay tetrahedralization.
    /// The ghost vertex is treated like a vertex at the center
    /// where all ghost tets have their orientations inverted.
    /// Before being calculated, this starts out as NAN.
    center: Pt3,
    /// Minimum distance between two nonadjacent simplexes of a tet.
    tolerance: f64,
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
    pub fn new(vertices: [(Pt3, V); 4], distance_tolerance: f64, default_tet: fn() -> T) -> Self
    where
        V: Clone,
        T: Clone,
    {
        let mut tets = Self {
            vertices: IdMap::default(),
            ghost_flags: Cell::new(VertexFlags::empty()),
            edges: IdMap::default(),
            tets: IdMap::default(),
            tri_map: FnvHashMap::default(),
            track_tri_map: false,
            center: Pt1::new(f64::NAN).xxx(),
            tolerance: distance_tolerance,
            default_tet,
        };

        tets.with_ids(
            [
                (VertexId(0), vertices[0].0, vertices[0].1.clone()),
                (VertexId(1), vertices[1].0, vertices[1].1.clone()),
                (VertexId(2), vertices[2].0, vertices[2].1.clone()),
                (VertexId(3), vertices[3].0, vertices[3].1.clone()),
            ],
            default_tet,
        );

        tets.center = tets.tet_centroid(TetId(0)).unwrap();

        tets
    }

    fn with_ids(&mut self, mut vertices: [(VertexId, Pt3, V); 4], default_tet: fn() -> T)
    where
        V: Clone,
        T: Clone,
    {
        // Make sure the tet is oriented positive
        vertices.sort_by_key(|(id, _, _)| *id);
        let swap = !sim::orient_3d(&vertices, |l, i| l[i].1.coords, 0, 1, 2, 3);

        let tw = |id, edge| TetWalker::new(TetId(id), edge);
        self.vertices
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

        let edges = [
            [vi[0], vi[1]],
            [vi[0], vi[2]],
            [vi[0], vi[3]],
            [vi[1], vi[2]],
            [vi[1], vi[3]],
            [vi[2], vi[3]],
            [vi[0], vi[4]],
            [vi[1], vi[4]],
            [vi[2], vi[4]],
            [vi[3], vi[4]],
        ];

        let ei = self
            .edges
            .extend_values(edges.iter().map(|vs| UndirEdge::new(*vs)))
            .collect::<Vec<_>>();

        let vertices_arr = [
            [vi[0], vi[1], vi[2], vi[3]],
            [vi[1], vi[2], vi[3], vi[4]],
            [vi[0], vi[3], vi[2], vi[4]],
            [vi[3], vi[0], vi[1], vi[4]],
            [vi[2], vi[1], vi[0], vi[4]],
        ];

        let edges_arr = [
            [ei[5], ei[0], ei[3], ei[2], ei[4], ei[1]],
            [ei[9], ei[3], ei[5], ei[7], ei[8], ei[4]],
            [ei[8], ei[2], ei[5], ei[6], ei[9], ei[1]],
            [ei[7], ei[2], ei[0], ei[9], ei[6], ei[4]],
            [ei[6], ei[3], ei[0], ei[8], ei[7], ei[1]],
        ];

        let opps_arr = [
            [tw(1, 7), tw(2, 7), tw(3, 7), tw(4, 7)],
            [tw(2, 8), tw(3, 5), tw(4, 2), tw(0, 4)],
            [tw(1, 8), tw(4, 5), tw(3, 2), tw(0, 5)],
            [tw(4, 8), tw(1, 5), tw(2, 2), tw(0, 6)],
            [tw(3, 8), tw(2, 5), tw(1, 2), tw(0, 7)],
        ];

        self.tets
            .extend_values(
                vertices_arr
                    .iter()
                    .zip(opps_arr.iter())
                    .zip(edges_arr.iter())
                    .map(|((vertices, opp_tets), edges)| {
                        Tet::new(*vertices, *opp_tets, *edges, default_tet())
                    }),
            )
            .for_each(|_| {});
    }

    /// Sets whether to track the triangle map.
    fn set_track_tri_map(&mut self, flag: bool) {
        if flag {
            self.track_tri_map = true;
            // Start now.
            self.tri_map.clear();
            for id in self.tets.keys() {
                for i in 0..4 {
                    let walker = TetWalker::new(id, i);
                    self.tri_map
                        .insert(TriVertices::new(walker.tri(self)), walker);
                }
            }
        } else {
            self.track_tri_map = false;
        }
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

    /// Gets a mutable reference to the flags of the vertex.
    /// If the ghost vertex is passed, this gets the ghost flags.
    fn vertex_flags_mut(&mut self, vertex: VertexId) -> &mut Cell<VertexFlags> {
        if vertex == Self::GHOST {
            &mut self.ghost_flags
        } else {
            &mut self.vertices[vertex].flags
        }
    }

    /// Gets the number of vertices in the tet mesh, not including the ghost vertex.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets the number of tets in the tet mesh, including ghost tets.
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

    /// Gets 6 times the volume of the tet. Returns None if one of its vertices is the ghost vertex.
    pub fn tet_volume_x6(&self, tet: TetId) -> Option<f64> {
        let vs = self[tet].vertices();
        if vs.contains(&Self::GHOST) {
            None
        } else {
            let pts = [
                self.vertices[vs[0]].position(),
                self.vertices[vs[1]].position(),
                self.vertices[vs[2]].position(),
                self.vertices[vs[3]].position(),
            ];

            Some(
                -(pts[1] - pts[0])
                    .cross(&(pts[2] - pts[0]))
                    .dot(&(pts[3] - pts[0])),
            )
        }
    }

    /// Gets the volume of the tet. Returns None if one of its vertices is the ghost vertex.
    pub fn tet_volume(&self, tet: TetId) -> Option<f64> {
        self.tet_volume_x6(tet).map(|v| v / 6.0)
    }

    /// Gets the centroid of the tet. Returns None if one of its vertices is the ghost vertex.
    pub fn tet_centroid(&self, tet: TetId) -> Option<Pt3> {
        let vs = self[tet].vertices();
        if vs.contains(&Self::GHOST) {
            None
        } else {
            let pts = [
                self.vertices[vs[0]].position(),
                self.vertices[vs[1]].position(),
                self.vertices[vs[2]].position(),
                self.vertices[vs[3]].position(),
            ];

            Some(((pts[0].coords + pts[1].coords + pts[2].coords + pts[3].coords) / 4.0).into())
        }
    }

    /// Gets a tet walker that starts at the given vertex.
    /// The vertex must be a solid vertex.
    pub fn walker_from_vertex(&self, vertex: VertexId) -> TetWalker {
        self.vertices[vertex].tet
    }

    /// Gets a tet walker that starts at the given edge, if one exists.
    /// The first vertex must be a solid vertex.
    pub fn walker_from_edge(&self, edge: [VertexId; 2]) -> Option<TetWalker> {
        self.vertex_target_walkers(edge[0])
            .find(|w| w.second(self) == edge[1])
    }

    /// Gets a tet walker that starts at the given triangle, if one exists.
    /// The first vertex must be a solid vertex.
    pub fn walker_from_tri(&self, tri: [VertexId; 3]) -> Option<TetWalker> {
        if self.track_tri_map {
            self.tri_map
                .get(&TriVertices::new(tri))
                .map(|w| w.to_edge_assuming_tri(self, [tri[0], tri[1]]))
        } else {
            self.walkers_from_vertex(tri[0])
                .flat_map(|w| {
                    if w.tri(self) == tri {
                        Some(w)
                    } else if w.to_nve().tri(self) == tri {
                        Some(w.to_nve())
                    } else if w.to_pve().tri(self) == tri {
                        Some(w.to_pve())
                    } else {
                        None
                    }
                })
                .next()
        }
    }

    /// Gets a tet walker at the given tetrahedron, if one exists.
    /// The first vertex must be a solid vertex.
    pub fn walker_from_tet(&self, tet: [VertexId; 4]) -> Option<TetWalker> {
        self.walker_from_tri([tet[0], tet[1], tet[2]]).filter(|w| w.fourth(self) == tet[3])
    }

    /// Gets a canonicalized tet walker that starts at the given tet.
    /// The tet must exist
    pub fn walker_from_tet_id(&self, tet: TetId) -> TetWalker {
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
    unsafe fn walkers_from_vertex_opt<'a>(
        &'a self,
        vertex: VertexId,
    ) -> WalkersFromVertexOpt<'a, V, T> {
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
        self.walkers_from_vertex_opt(vertex)
            .map(|walker| walker.id())
    }

    /// Gets the tet walkers for the edges from a vertex.
    pub fn vertex_target_walkers<'a>(&'a self, vertex: VertexId) -> VertexTargetWalkers<'a, V, T> {
        VertexTargetWalkers {
            mesh: self,
            visited: FnvHashSet::default(),
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the tet walkers for the edges from a vertex.
    ///
    /// # Safety
    /// The returned iterator needs to drop before calling this function again.
    unsafe fn vertex_target_walkers_opt<'a>(
        &'a self,
        vertex: VertexId,
    ) -> VertexTargetWalkersOpt<'a, V, T> {
        VertexTargetWalkersOpt {
            mesh: self,
            visited: vec![],
            to_search: vec![self.walker_from_vertex(vertex)],
        }
    }

    /// Gets the edge target vertices from a vertex.
    pub fn vertex_targets<'a>(&'a self, vertex: VertexId) -> VertexTargets<'a, V, T> {
        VertexTargets(self.vertex_target_walkers(vertex))
    }

    /// Gets the edge target vertices from a vertex.
    ///
    /// # Safety
    /// The returned iterator needs to drop before calling this function again.
    unsafe fn vertex_targets_opt<'a>(&'a self, vertex: VertexId) -> VertexTargetsOpt<'a, V, T> {
        VertexTargetsOpt(self.vertex_target_walkers_opt(vertex))
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
            if v1 == Self::GHOST {
                0.0
            } else {
                f64::INFINITY
            }
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

    /// The indexing function to use for orient_3d.
    fn idx_ori(&self, vertex: VertexId) -> Vec3 {
        if vertex == Self::GHOST {
            self.center.coords
        } else {
            self[vertex].position().coords
        }
    }

    fn idx(&self, vertex: VertexId) -> Vec3 {
        self[vertex].position().coords
    }

    /// Gets whether these points are oriented positive.
    /// Uses simulation of simplicity to avoid ties.
    /// Any vertex is allowed to be the ghost vertex and this should just work™.
    pub fn orient_3d(&self, v0: VertexId, v1: VertexId, v2: VertexId, v3: VertexId) -> bool {
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
    pub fn in_sphere_sim(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
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
    pub fn tri_intersects_edge(
        &self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
        v4: VertexId,
    ) -> bool {
        let positive = self.orient_3d(v2, v1, v0, v3);
        self.orient_3d(v0, v1, v2, v4) == positive &&
        self.orient_3d(v0, v1, v3, v4) == positive &&
        self.orient_3d(v1, v2, v3, v4) == positive &&
        self.orient_3d(v2, v0, v3, v4) == positive
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
        self.orient_3d(v0, v1, v2, v4) &&
        self.orient_3d(v3, v2, v1, v4) &&
        self.orient_3d(v2, v3, v0, v4) &&
        self.orient_3d(v1, v0, v3, v4)
    }

    /// Flips in a vertex using one big mega flip.
    /// The boundary must be manifold.
    /// For now, no vertex can be completely enclosed.
    /// This does not check that negative-volume tets won't be created. Be careful.
    pub fn flip_in_vertex_unchecked(
        &mut self,
        vertex: VertexId,
        boundary: &mut [TetWalker],
        enclosure: &[TetId],
    ) where
        T: Clone,
    {
        if self.track_tri_map {
            for id in enclosure {
                for i in 0..4 {
                    self.tri_map
                        .remove(&TriVertices::new(TetWalker::new(*id, i).tri(self)));
                }
            }
        }

        // Delete tets that don't even have a triangle on the boundary
        // Mark boundary
        for walker in &*boundary {
            self.tets[walker.id()].set_flags(TetFlags::BOUNDARY);
            self.edges[walker.undir_edge(self)].set_flags(UndirEdgeFlags::BOUNDARY);
            self.edges[walker.to_nfe().undir_edge(self)].set_flags(UndirEdgeFlags::BOUNDARY);
            self.edges[walker.to_pfe().undir_edge(self)].set_flags(UndirEdgeFlags::BOUNDARY);
        }
        // Removed unmarked enclosure
        for id in enclosure {
            for edge in &self.tets[*id].edges {
                if self
                    .edges
                    .get(*edge)
                    .map(|e| !e.flags.get().intersects(UndirEdgeFlags::BOUNDARY))
                    .unwrap_or(false)
                {
                    self.edges.remove(*edge);
                }
            }

            if !self.tets[*id].flags.get().intersects(TetFlags::BOUNDARY) {
                self.tets.remove(*id);
            }
        }

        // Check for tets that need to be duplicated and make necessary clones
        for walker in boundary.iter_mut() {
            if self.tets[walker.id()]
                .flags
                .get()
                .intersects(TetFlags::NEEDS_DUPLICATION)
            {
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
        // Adjacency map is now stored in edges.

        for walker in &*boundary {
            // Also unmark boundary and tet duplication flags
            self.tets[walker.id()].clear_flags(TetFlags::BOUNDARY | TetFlags::NEEDS_DUPLICATION);
            self.edges[walker.undir_edge(self)].clear_flags(UndirEdgeFlags::BOUNDARY);
            self.edges[walker.to_nfe().undir_edge(self)].clear_flags(UndirEdgeFlags::BOUNDARY);
            self.edges[walker.to_pfe().undir_edge(self)].clear_flags(UndirEdgeFlags::BOUNDARY);

            // Don't forget the new vertex!
            self.tets[walker.id()].vertices[walker.edge as usize % 4] = vertex;
            if self[vertex].tet.id() == TetId::invalid() {
                self.vertices[vertex].tet = walker.to_perm(Permutation::_3210);
            }
            // Other vertices, just in case their references are now dangling
            let tri = walker.tri(self);
            self.vertices.get_mut(tri[0]).map(|v| v.tet = *walker);
            self.vertices
                .get_mut(tri[1])
                .map(|v| v.tet = walker.to_nfe());
            self.vertices
                .get_mut(tri[2])
                .map(|v| v.tet = walker.to_pfe());

            // Edge adjacency and internal links
            for walker in &[*walker, walker.to_nfe(), walker.to_pfe()] {
                let edge = walker.undir_edge(self);
                let boundary_walker = &mut self.edges[edge].boundary_walker;
                if boundary_walker.id().is_valid() {
                    let walker = walker.to_twin_edge();
                    let adj = boundary_walker.to_twin_edge();

                    // Calling `walker.to_adj(&self)` should give `adj`
                    self.tets[walker.id()].opp_tets[walker.edge as usize % 4] =
                        match walker.edge / 4 {
                            0 => adj,
                            1 => adj.to_nfe(),
                            2 => adj.to_pfe(),
                            _ => unreachable!(),
                        };

                    // Symmetrically
                    self.tets[adj.id()].opp_tets[adj.edge as usize % 4] = match adj.edge / 4 {
                        0 => walker,
                        1 => walker.to_nfe(),
                        2 => walker.to_pfe(),
                        _ => unreachable!(),
                    };
                    boundary_walker.tet = TetId::invalid();
                } else {
                    *boundary_walker = *walker;
                }
            }
        }

        // Add edges.
        for walker in &*boundary {
            for perm in &[Permutation::_0312, Permutation::_1320, Permutation::_2301] {
                let walker = walker.to_perm(*perm);
                let source = walker.first(self);
                if !self
                    .vertex_flags(source)
                    .get()
                    .intersects(VertexFlags::EDGE_ADDED)
                {
                    *self.vertex_flags_mut(source).get_mut() |= VertexFlags::EDGE_ADDED;

                    let edge = self.edges.insert(UndirEdge::new([source, vertex]));
                    let mut iter = walker;
                    while {
                        self.tets[iter.id()].edges[iter.undir_edge_index()] = edge;
                        iter = iter.to_twin_edge().to_adj(self);
                        iter != walker
                    } {}
                }
            }
        }

        for walker in &*boundary {
            *self.vertex_flags_mut(walker.first(self)).get_mut() &= !VertexFlags::EDGE_ADDED;
            *self.vertex_flags_mut(walker.to_nfe().first(self)).get_mut() &=
                !VertexFlags::EDGE_ADDED;
            *self.vertex_flags_mut(walker.to_pfe().first(self)).get_mut() &=
                !VertexFlags::EDGE_ADDED;
        }

        if self.track_tri_map {
            for walker in &*boundary {
                for i in 0..4 {
                    let walker = walker.to_edge(i);
                    self.tri_map
                        .insert(TriVertices::new(walker.tri(self)), walker);
                }
            }
        }
    }

    /// Create a Delaunay tetrahedralization from vertices. Panics if there are less than 4 vertices
    /// because it takes 4 vertices to make a tet.
    /// `VertexId(i)` is the index to the `i`th vertex returned from the iterator.
    pub fn delaunay_from_vertices<I: IntoIterator<Item = (Pt3, V)>>(
        vertices: I,
        distance_tolerance: f64,
        default_tet: fn() -> T,
    ) -> Self
    where
        V: Clone,
        T: Clone,
    {
        let mut mesh = Self {
            vertices: IdMap::default(),
            ghost_flags: Cell::new(VertexFlags::empty()),
            edges: IdMap::default(),
            tets: IdMap::default(),
            tri_map: FnvHashMap::default(),
            track_tri_map: false,
            center: Pt1::new(f64::NAN).xxx(),
            tolerance: distance_tolerance,
            default_tet,
        };
        let mut vertices = util::hilbert_sort(vertices);

        // Pre-add the vertices to keep the id map dense
        mesh.vertices
            .extend_values((0..vertices.len()).map(|_| {
                Vertex::new(
                    TetWalker::new(TetId::invalid(), 0),
                    Pt3::origin(),
                    vertices[0].2.clone(),
                )
            }))
            .for_each(|_| {});

        // Set the first 4 points to not be coplanar.
        // Assume the first 2 points are not the same because of the distance tolerance.
        let p0 = vertices[0].1;
        let p1 = vertices[1].1;
        let v2i = vertices.iter().enumerate().skip(2).map(|(i, (_, p2, _))| {
            let areax2 = (p1 - p0).cross(&(p2 - p0)).norm();
            let max_len = (p1 - p0).norm().max((p2 - p1).norm()).max((p0 - p2).norm());
            (i, areax2 / max_len)
        }).filter(|(_, r)| *r >= distance_tolerance).max_by_key(|(_, r)| FloatOrd(*r))
            .expect("All points are collinear within the tolerance.").0;
        vertices.swap(2, v2i);
        let p2 = vertices[2].1;

        let v3i = vertices.iter().skip(3).position(|(_, p3, _)| {
            let volx6 = -(p1 - p0).cross(&(p2 - p0)).dot(&(p3 - p0));
            let max_area = [[p0, p1, p0, p2], [*p3, p2, *p3, p1], [p2, *p3, p2, p0], [p1, p0, p1, *p3],
                [p0, p1, p2, *p3], [p0, p2, p1, *p3], [p0, *p3, p1, p2]].iter().map(|[pa, pb, pc, pd]|
                    (pb - pa).cross(&(pd - pc)).norm())
                    .max_by_key(|area| FloatOrd(*area)).unwrap();
            volx6.abs() / max_area >= distance_tolerance
        }).expect("All points are coplanar within the tolerance.");
        vertices.swap(3, v3i + 3);

        let mut drain = vertices.drain(0..4);
        let v0 = drain.next().unwrap();
        let v1 = drain.next().unwrap();
        let v2 = drain.next().unwrap();
        let v3 = drain.next().unwrap();
        std::mem::drop(drain);
        let v3_id = v3.0;

        mesh.with_ids([v0, v1, v2, v3], default_tet);
        // include v3's id to avoid indexing out of bounds
        let ids = iter::once(v3_id)
            .chain(vertices.iter().map(|v| v.0))
            .collect::<Vec<_>>();

        for (i, (id, position, value)) in vertices.into_iter().enumerate() {
            mesh.add_vertex_delaunay_internal(Some(id), position, value, ids[i]);
        }

        // Obtain center of mass. The axis-aligned bounding box's center is not always in the tetrahedralization.
        let mut weight_sum = 0.0;
        let mut center_sum = Vec3::from_element(0.0);
        for tet in mesh.tets.keys() {
            let weight = mesh.tet_volume_x6(tet).unwrap_or(0.0);
            let center = mesh.tet_centroid(tet).unwrap_or(Pt3::origin());
            weight_sum += weight;
            center_sum += center.coords * weight;
        }

        mesh.center = (center_sum / weight_sum).into();

        mesh
    }

    /// Add a vertex to keep the Delaunay property, assuming this is a Delaunay tetrahedralization.
    /// Returns the new vertex id.
    pub fn add_vertex_delaunay(&mut self, position: Pt3, value: V) -> VertexId
    where
        T: Clone,
    {
        self.add_vertex_delaunay_internal(
            None,
            position,
            value,
            self.vertices.keys().next().unwrap(),
        )
    }

    pub fn add_vertex_delaunay_internal(
        &mut self,
        id: Option<VertexId>,
        position: Pt3,
        value: V,
        search_start: VertexId,
    ) -> VertexId
    where
        T: Clone,
    {
        let mut closest = search_start;

        let id = if let Some(id) = id {
            self.vertices.insert_with_key(
                id,
                Vertex::new(TetWalker::new(TetId::invalid(), 0), position, value),
            );
            id
        } else {
            self.vertices.insert(Vertex::new(
                TetWalker::new(TetId::invalid(), 0),
                position,
                value,
            ))
        };

        let lucky_tet = self[closest].tet.id();
        let not_delaunay = if {
            let vs = self[lucky_tet].vertices;
            self.in_sphere_sim(vs[0], vs[1], vs[2], vs[3], id)
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
                    self.in_sphere_sim(vs[0], vs[1], vs[2], vs[3], id)
                })
            }
            .unwrap_or_else(|| {
                // Rare case because the distance comparisons are not robust
                // Safe because the previous call's iterator was dropped.
                let mut checked =
                    unsafe { self.vertex_tets_opt(closest).collect::<FnvHashSet<_>>() }; // already checked
                let mut tets = checked
                    .iter()
                    .flat_map(|tet| self.adjacent_tets(*tet).to_vec())
                    .collect::<FnvHashSet<_>>()
                    .into_iter()
                    .collect::<VecDeque<_>>();

                while let Some(tet) = tets.pop_front() {
                    if checked.insert(tet) {
                        let vs = self[tet].vertices();
                        if self.in_sphere_sim(vs[0], vs[1], vs[2], vs[3], id) {
                            return tet;
                        }
                    }
                }

                panic!(
                    "Expected some tet to stop being Delaunay when inserting vertex {} at {:?}",
                    id,
                    self[id].position()
                )
            })
        };

        // Find all non-Delaunay tets
        let (mut boundary, enclosed) =
            self.walker_from_tet_id(not_delaunay)
                .point_cavity(self, true, id);

        // Add the new vertex
        self.flip_in_vertex_unchecked(id, &mut boundary, &enclosed);

        id
    }

    /// Improves the quality of the tet mesh with repeated flips.
    pub fn flip_quality(&mut self) where V: Clone, T: Clone {
        let mut to_flip = self.edges.values().map(|edge| EdgeOrTri::Edge(sorted_2(edge.vertices)))
            .chain(self.tets.keys().flat_map(|id| vec![
                TetWalker::new(id, 0),
                TetWalker::new(id, 1),
                TetWalker::new(id, 2),
                TetWalker::new(id, 3),
            ]).map(|w| w.tri(self))
            .filter(|tri| even_sorted_3(*tri) == sorted_3(*tri))
            .map(EdgeOrTri::Tri))
            .collect::<VecDeque<_>>();
        let mut flip_set = to_flip.iter().copied().collect::<FnvHashSet<_>>();
        
        while let Some(simplex) = to_flip.pop_front() {
            flip_set.remove(&simplex);

            match simplex {
                EdgeOrTri::Edge(edge) => {
                    if let Some(walker) = self.walker_from_edge(edge) {
                        if let Some(walkers) = walker.flip32_quality(self) {
                            for edge_or_tri in &[
                                EdgeOrTri::Tri(sorted_3(walkers[0].tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[1].tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[0].to_nve().tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[1].to_nve().tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[0].to_pve().tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[1].to_pve().tri(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[0].edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[1].edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[0].to_nve().edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[1].to_nve().edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[0].to_pve().edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[1].to_pve().edge(self))),
                            ] {
                                if flip_set.insert(*edge_or_tri) {
                                    to_flip.push_back(*edge_or_tri);
                                }
                            }
                        }
                    }
                }

                EdgeOrTri::Tri(tri) => {
                    // tri should be sorted.
                    if let Some(walker) = self.tri_map.get(&TriVertices(tri)) {
                        if let Some(walkers) = walker.flip23_quality(self) {
                            for edge_or_tri in &[
                                EdgeOrTri::Tri(sorted_3(walkers[0].tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[1].tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[2].tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[0].to_twin_edge().tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[1].to_twin_edge().tri(self))),
                                EdgeOrTri::Tri(sorted_3(walkers[2].to_twin_edge().tri(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[0].edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[1].edge(self))),
                                EdgeOrTri::Edge(sorted_2(walkers[2].edge(self))),
                            ] {
                                if flip_set.insert(*edge_or_tri) {
                                    to_flip.push_back(*edge_or_tri);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Removes silver tets in the mesh.
    pub fn remove_slivers(&mut self) where V: Clone, T: Clone {
        let slivers = self.tets.iter().filter(|(id, tet)| {
            let vs = tet.vertices;
            !self.non_sliver(vs[0], vs[1], vs[2], vs[3])
        }).map(|(_, tet)| tet.vertices).collect::<Vec<_>>();

        println!("Sliver to tet ratio: {}/{}", slivers.len(), self.num_tets());

        let mut num_slivers = slivers.len();
        for sliver in slivers {
            for depth in &[1, 2, 4, 8] {
                if let Some(walker) = self.walker_from_tet(sliver) {
                    match walker.sliver_kind(self) {
                        SliverKind::EdgeEdge(walker) => {
                            let other_tet = walker.to_opp_edge().tet(self);
                            if walker.progress_by_removing_edge(
                                self, |_, _| true, |_, _| true, *depth, false, &mut vec![], &mut vec![]
                            ).is_some() {
                                num_slivers -= 1;
                                break;
                            }

                            // Walker got invalidated
                            let walker = self.walker_from_tet(other_tet).unwrap();
                            if walker.progress_by_removing_edge(
                                self, |_, _| true, |_, _| true, *depth, false, &mut vec![], &mut vec![]
                            ).is_some() {
                                num_slivers -= 1;
                                break;
                            }
                        }

                        SliverKind::TriVertex(walker) => {
                            if walker.progress_by_removing_tri(
                                self, |_, _| true, |_, _| true, *depth, false, &mut vec![], &mut vec![]
                            ).is_some() {
                                num_slivers -= 1;
                                break;
                            }
                        }
                    }
                }
            }
        }

        println!("New sliver to tet ratio: {}/{}", num_slivers, self.num_tets());
        self.export_debug_obj("ignore/less_slivers.obj").unwrap();
    }

    /// Attempts to recover an edge
    /// and locks it to prevent it from flipping outside of unchecked flip functions.
    /// Returns whether it succeeded.
    ///
    /// Expects both vertices to be solid vertices.
    pub fn recover_and_lock_edge(&mut self, edge: [VertexId; 2], depth: usize) -> bool where V: Clone, T: Clone {
        // Don't store walkers because they get invalidated by flips.
        let mut intersecting: Vec<[VertexId; 3]> = vec![];

        loop {
            if let Some(tri) = intersecting.last().copied() {
                // Look for an intersecting triangle from this one, which couldn't be removed
                let walker = self.walker_from_tri(tri).unwrap().to_adj_ae(self);
                if walker.fourth(self) == edge[1] {
                    return false;
                }

                let mut reached = false;
                for walker in &[walker.to_twin_edge(), walker.to_opp_edge(), walker.to_perm(Permutation::_3210)] {
                    let tri = walker.tri(self);
                    if self.tri_intersects_edge(tri[0], tri[1], tri[2], edge[0], edge[1]) {
                        intersecting.push(tri);
                        reached = true;
                        break;
                    }
                }
                assert!(reached);
            } else {
                // Look for an intersecting triangle from the vertex.
                if let Some(tri) = self.walkers_from_vertex(edge[0])
                    .map(|w| w.opp_tri(self))
                    .find(|tri| {
                        !tri.contains(&edge[1]) &&
                        self.tri_intersects_edge(tri[0], tri[1], tri[2], edge[0], edge[1])
                    })
                {
                    intersecting.push(tri);
                } else {
                    // No triangle intersects the edge, so the edge must exist now.
                    let walker = self.walker_from_edge(edge).unwrap();
                    self.edges[walker.undir_edge(self)].set_flags(UndirEdgeFlags::LOCKED);
                    return true;
                }
            };

            // Remove as many intersecting triangles as possible until blocked.
            while let Some(tri) = intersecting.pop() {
                if let Some(walker) = self.walker_from_tri(tri) {
                    if let Some(goal_index) = walker.progress_by_removing_tri(
                        self,
                        |_, _| true,
                        |mesh, tri| 
                            tri.contains(&edge[0]) || tri.contains(&edge[1]) ||
                            !mesh.tri_intersects_edge(tri[0], tri[1], tri[2], edge[0], edge[1]),
                        depth,
                        true,
                        &mut vec![],
                        &mut vec![],
                    ) {
                        assert_eq!(goal_index, 0);
                    } else {
                        // get blokt
                        intersecting.push(tri);
                        break;
                    }
                }
            }
        }
    }

    /// Recovers specific edges in a PLC, adding Steiner points as necessary in both the tet mesh and the PLC,
    /// and locks them to prevent them from flipping outside of unchecked flip functions.
    pub fn recover_and_lock_specific_edges<I: IntoIterator<Item = [VertexId; 2]>>(&mut self, edges: I) -> usize where V: Clone, T: Clone {
        let mut to_recover = edges.into_iter().collect::<Vec<_>>();
        let initial = to_recover.clone();
        let mut curr_to_recover = vec![];
        let mut num_steiner = 0;
        while !to_recover.is_empty() {
            for depth in &[1, 2, 4, 8] {
                curr_to_recover.append(&mut to_recover);
                for edge in curr_to_recover.drain(..) {
                    if !self.recover_and_lock_edge(edge, *depth) {
                        to_recover.push(edge);
                    }
                }
            }

            curr_to_recover.append(&mut to_recover);
            for edge in curr_to_recover.drain(..) {
                if self.walker_from_edge(edge).is_none() {
                    // Look for an intersecting triangle from the vertex.
                    let start = self.walkers_from_vertex(edge[0])
                        .find(|walker| {
                            let opp = walker.opp_tri(self);
                            self.tri_intersects_edge(opp[0], opp[1], opp[2], edge[0], edge[1])
                        }).unwrap().to_perm(Permutation::_3210);

                    // Look for all tets that intersect the edge.
                    let mut end = false;
                    let search = iter::successors(Some(start), |walker| {
                        let adj = walker.to_adj_ae(self);
                        if end {
                            None
                        } else if adj.fourth(self) == edge[1] {
                            end = true;
                            Some(adj)
                        } else {
                            Some(*[adj.to_twin_edge(), adj.to_opp_edge(), adj.to_perm(Permutation::_3210)].iter().find(|walker| {
                                let tri = walker.tri(self);
                                self.tri_intersects_edge(tri[0], tri[1], tri[2], edge[0], edge[1])
                            }).unwrap())
                        }
                    }).collect::<VecDeque<_>>();

                    // Add a Steiner point at the midpoint.
                    let mut search_iter = search.iter();
                    search_iter.next_back(); // Don't include the last element

                    //println!();
                    let ep = [self[edge[0]].position(), self[edge[1]].position()];
                    let mid = search_iter.map(|tet| {
                        // Line-plane intersection
                        let tri = tet.tri(self);
                        let tp = [self[tri[0]].position(), self[tri[1]].position(), self[tri[2]].position()];

                        let numer = -robust_geo::orient_3d(ep[0].coords, tp[1].coords, tp[2].coords, tp[0].coords);
                        let denom = robust_geo::orient_3d(ep[1] - ep[0] + tp[0].coords, tp[1].coords, tp[2].coords, tp[0].coords);
                        numer / denom
                    })
                    //.inspect(|fac| {
                    //    if initial[0][0].0 == 1871 && initial[0][1].0 == 2014 {
                    //        println!("Len factor: {}", (ep[1] - ep[0]).norm() * fac);
                    //    }
                    //})
                    //.inspect(|(n, d)| {
                    //    println!("{}/{}", n, d);
                    //    let len = (ep[1] - ep[0]).norm();
                    //    if len * n / d < self.tolerance || len * n / d > len - self.tolerance {
                    //        self.export_debug_obj("ignore/debug_obj.obj").unwrap();
                    //        panic!("Factor too close to 0 or 1, inserting in edge {}-{}, position {:?} to {:?}", edge[0], edge[1], ep[0], ep[1]);
                    //    }
                    //})
                    .filter(|fac| {
                        let len = (ep[1] - ep[0]).norm();
                        len * *fac >= self.tolerance && len * *fac <= len - self.tolerance
                        //*fac >= 0.0 && *fac <= 1.0
                    })
                    .min_by_key(|fac| FloatOrd((fac - 0.5).abs()))
                    .map_or(((ep[0].coords + ep[1].coords) / 2.0).into(),
                        |fac| ep[0] + (ep[1] - ep[0]) * fac);
                    //let mid = Pt3::from((ep[0].coords + ep[1].coords) / 2.0);
                    //if initial[0][0].0 == 1871 && initial[0][1].0 == 2014 {
                    //    println!();
                    //}
                    assert!((mid - ep[0]).norm() >= self.tolerance && (mid - ep[1]).norm() >= self.tolerance,
                        "Steiner point overload when trying to recover edge {}-{}", initial[0][0], initial[0][1]);

                    let steiner = self.vertices.insert(
                        Vertex::new(TetWalker::new(TetId::invalid(), 0), mid, self[edge[0]].value().clone()));
                    
                    let mut search = search.into_iter().map(TetWalker::id).collect::<VecDeque<_>>();

                    let mut visited = FnvHashSet::default();
                    while let Some(tet) = search.pop_front() {
                        if visited.insert(tet) {
                            let vs = self[tet].vertices();
                            let walker = self.walker_from_tet_id(tet);
                            if self.tet_intersects_vertex(vs[0], vs[1], vs[2], vs[3], steiner) {
                                let (mut boundary, enclosure) = walker.point_cavity(self, false, steiner);
                                self.flip_in_vertex_unchecked(steiner, &mut boundary, &enclosure);
                                to_recover.push([edge[0], steiner]);
                                to_recover.push([steiner, edge[1]]);
                                num_steiner += 1;
                                break;
                            } else {
                                search.push_back(walker.to_adj_ae(self).id());
                                search.push_back(walker.to_twin_edge().to_adj_ae(self).id());
                                search.push_back(walker.to_opp_edge().to_adj_ae(self).id());
                                search.push_back(walker.to_perm(Permutation::_3210).to_adj_ae(self).id());
                            }
                        }

                        assert!(!search.is_empty(), "Adding Steiner point {} in edge {}-{} failed", steiner, edge[0], edge[1]);
                    }
                }
            }
        }

        num_steiner
    }

    /// Recovers edges in a PLC, adding Steiner points as necessary in both the tet mesh and the PLC,
    /// and locks them to prevent them from flipping outside of unchecked flip functions.
    pub fn recover_and_lock_edges<E, F>(&mut self, plc: &mut Plc<V, E, F>) where V: Clone, T: Clone {
        // Needed for the edge and tri removal functions.
        self.set_track_tri_map(true);

        self.remove_slivers();

        let mut plc_edges = plc
            .edges()
            .map(|(_, edge)| sorted_2(edge.vertices()))
            .collect::<FnvHashSet<_>>();

        // Lock already recovered edges
        for (_, edge) in self.edges.iter() {
            if plc_edges.remove(&edge.vertices) {
                // No need to remove the reverse because `edge.vertices` is sorted.
                edge.set_flags(UndirEdgeFlags::LOCKED);
            }
        }

        let to_recover = plc_edges.into_iter().collect::<Vec<_>>();
        let mut total_steiner = 0;
        for edge in to_recover {
            let num_steiner = self.recover_and_lock_specific_edges(iter::once(edge));
            if num_steiner > 0 {
                println!("Added {} Steiner points to recover {}-{}", num_steiner, edge[0], edge[1]);
            }
            total_steiner += num_steiner;
        }
        println!("Added {} Steiner points total.", total_steiner);
    }

    #[cfg(feature = "obj")]
    pub fn export_debug_obj<P: AsRef<Path>>(&self, path: P) -> Result<(), obj::ObjError> {
        let solid_tets = || {
            self.tets
                .values()
                .filter(|t| !t.vertices().contains(&Self::GHOST))
                .enumerate()
        };

        let obj = obj::ObjData {
            position: solid_tets()
                .flat_map(|(_, t)| t.vertices().to_vec())
                .map(|v| {
                    let pos = self[v].position;
                    [pos.x as f32, pos.y as f32, pos.z as f32]
                })
                .collect(),
            texture: vec![],
            normal: vec![],

            objects: vec![obj::Object {
                name: "Tet Mesh".to_owned(),
                groups: vec![obj::Group {
                    name: "Tet Group".to_owned(),
                    index: 0,
                    material: None,
                    polys: solid_tets()
                        .flat_map(|(i, _)| {
                            vec![
                                [4 * i + 0, 4 * i + 1, 4 * i + 2],
                                [4 * i + 3, 4 * i + 2, 4 * i + 1],
                                [4 * i + 2, 4 * i + 3, 4 * i + 0],
                                [4 * i + 1, 4 * i + 0, 4 * i + 3],
                            ]
                        })
                        .map(|[a, b, c]| {
                            obj::SimplePolygon(vec![
                                obj::IndexTuple(a, None, None),
                                obj::IndexTuple(b, None, None),
                                obj::IndexTuple(c, None, None),
                            ])
                        })
                        .collect::<Vec<_>>(),
                }],
            }],

            material_libs: vec![],
        };
        obj.save(path)
    }

    /// Assert all the invariants of this tet mesh.
    #[track_caller]
    fn assert_integrity(self: &TetMesh<V, T>, orient_test: bool) {
        // Vertex invariants
        for (id, vertex) in self.vertices.iter() {
            assert!(
                vertex.tet.id().is_valid(),
                "Vertex {} has an invalid tet walker",
                id
            );
            assert_eq!(
                vertex.tet.first(self),
                id,
                "Tet walker ({}, {}) of vertex {} does not start with the vertex.",
                vertex.tet.id(),
                vertex.tet.edge,
                id
            );
        }

        // Edge invariants
        for (id, edge) in self.edges.iter() {
            assert_ne!(
                edge.vertices[0], edge.vertices[1],
                "Edge {} has equal vertices.",
                id
            );
            assert!(
                edge.vertices[0] < edge.vertices[1],
                "Edge {} has unsorted vertices.",
                id
            );
        }
        let edge_set = self
            .edges
            .values()
            .map(|edge| {
                let mut vs = edge.vertices;
                vs.sort();
                vs
            })
            .collect::<FnvHashSet<_>>();
        assert_eq!(
            self.edges.len(),
            edge_set.len(),
            "There are duplicate edges."
        );

        let mut edge_ids = self.edges.keys().collect::<FnvHashSet<_>>();
        for tet in self.tets.values() {
            for id in &tet.edges {
                edge_ids.remove(id);
            }
        }
        if !edge_ids.is_empty() {
            panic!(
                "Edge {} does not belong to a tet.",
                edge_ids.iter().next().unwrap()
            );
        }

        if self.track_tri_map {
            for (tri, walker) in &self.tri_map {
                assert_eq!(
                    TriVertices::new(walker.tri(self)),
                    *tri,
                    "Triangle {}-{}-{} points to the wrong tet walker: ({}, {})",
                    tri.0[0],
                    tri.0[1],
                    tri.0[2],
                    walker.id(),
                    walker.edge
                );
            }
        }

        // Tet invariants
        for (id, tet) in self.tets.iter() {
            assert_eq!(
                tet.vertices().len(),
                tet.vertices().iter().collect::<FnvHashSet<_>>().len(),
                "Tet {} does not have unique vertices.",
                id
            );

            if orient_test {
                let vs = tet.vertices();
                assert!(
                    self.orient_3d(vs[0], vs[1], vs[2], vs[3]),
                    "Tet {} (vertices {} {} {} {}) is oriented negative.",
                    id,
                    vs[0],
                    vs[1],
                    vs[2],
                    vs[3]
                );
            }

            for i in 0..12 {
                let walker = TetWalker::new(id, i);
                let vi = walker.tet(self);
                let mut vs = vi;
                vs.swap(0, 1);
                vs[3] = tet.opp_tets[i as usize % 4].fourth(self);
                assert_eq!(walker.to_adj(self).tet(self), vs, 
                    "Tet walker ({}, {}) (vertices {} {} {} {}) does not point to the adjacent tet the right way",
                    walker.id(), walker.edge, vi[0], vi[1], vi[2], vi[3]);

                assert_eq!(walker.to_adj(self).to_adj(self).tet(self), vi,
                    "Adjacent tet walker ({}, {}) of tet walker ({}, {}) (vertices {} {} {} {}) does not point back to the tet walker.",
                    walker.to_adj(self).id(), walker.to_adj(self).edge, walker.id(), walker.edge, vi[0], vi[1], vi[2], vi[3]);

                let mut edge_vs = self.edges[walker.undir_edge(self)].vertices;
                vs[..2].sort();
                edge_vs.sort();
                assert_eq!(
                    vs[..2],
                    edge_vs,
                    "Tet walker ({}, {}) does not point to the right edge.",
                    walker.id(),
                    walker.edge
                );

                if i < 4 && self.track_tri_map {
                    assert!(
                        self.tri_map
                            .contains_key(&TriVertices::new(walker.tri(self))),
                        "Tet walker ({}, {}) has triangle missing from map",
                        walker.id(),
                        walker.edge
                    );
                }
            }
        }
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
                self.to_search
                    .push(walker.to_perm(Permutation::_1032).to_adj(self.mesh));
                self.to_search
                    .push(walker.to_perm(Permutation::_2013).to_adj(self.mesh));
                self.to_search
                    .push(walker.to_perm(Permutation::_3021).to_adj(self.mesh));
                return Some(walker);
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
                self.to_search
                    .push(walker.to_perm(Permutation::_1032).to_adj(self.mesh));
                self.to_search
                    .push(walker.to_perm(Permutation::_2013).to_adj(self.mesh));
                self.to_search
                    .push(walker.to_perm(Permutation::_3021).to_adj(self.mesh));
                return Some(walker);
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
pub struct VertexTargetWalkers<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: FnvHashSet<VertexId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for VertexTargetWalkers<'a, V, T> {
    type Item = TetWalker;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(walker) = self.to_search.pop() {
            let target = walker.second(self.mesh);

            if self.visited.insert(target) {
                // Search is different; search 1 edge radius at a time.
                let mut adj = walker.to_adj(self.mesh).to_nfe();
                let start = adj;
                while {
                    self.to_search.push(adj);
                    adj = adj.to_perm(Permutation::_3021).to_adj(self.mesh);
                    adj != start
                } {}

                return Some(walker);
            }
        }
        None
    }
}

#[derive(Clone, Debug)]
pub struct VertexTargetWalkersOpt<'a, V, T> {
    mesh: &'a TetMesh<V, T>,
    visited: Vec<VertexId>,
    to_search: Vec<TetWalker>,
}

impl<'a, V, T> Iterator for VertexTargetWalkersOpt<'a, V, T> {
    type Item = TetWalker;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(walker) = self.to_search.pop() {
            let target = walker.second(self.mesh);

            let flags = self.mesh.vertex_flags(target);
            if !flags.get().intersects(VertexFlags::VERTEX_TARGETS) {
                flags.set(flags.get() | VertexFlags::VERTEX_TARGETS);
                self.visited.push(target);

                // Search is different; search 1 edge radius at a time.
                let mut adj = walker.to_adj(self.mesh).to_nfe();
                let start = adj;
                while {
                    self.to_search.push(adj);
                    adj = adj.to_perm(Permutation::_3021).to_adj(self.mesh);
                    adj != start
                } {}

                return Some(walker);
            }
        }
        None
    }
}

impl<'a, V, T> Drop for VertexTargetWalkersOpt<'a, V, T> {
    fn drop(&mut self) {
        for target in &self.visited {
            let flags = self.mesh.vertex_flags(*target);
            flags.set(flags.get() & !VertexFlags::VERTEX_TARGETS);
        }
    }
}

pub struct VertexTargets<'a, V, T>(VertexTargetWalkers<'a, V, T>);

impl<'a, V, T> Iterator for VertexTargets<'a, V, T> {
    type Item = VertexId;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|w| w.second(self.0.mesh))
    }
}

struct VertexTargetsOpt<'a, V, T>(VertexTargetWalkersOpt<'a, V, T>);

impl<'a, V, T> Iterator for VertexTargetsOpt<'a, V, T> {
    type Item = VertexId;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|w| w.second(self.0.mesh))
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
        let mut tets = TetMesh::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
            ],
            0.0,
            || (),
        );
        tets.set_track_tri_map(true);
        tets
    }

    fn reverse_tets() -> TetMesh<(), ()> {
        let mut tets = TetMesh::new(
            [
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
            ],
            0.0,
            || (),
        );
        tets.set_track_tri_map(true);
        tets
    }

    #[track_caller]
    fn assert_tets<V, T, I: IntoIterator<Item = [u32; 4]>>(mesh: &TetMesh<V, T>, tets: I) {
        let result = mesh
            .tets
            .iter()
            .map(|(_, tet)| even_sort_4(tet.vertices()))
            .collect::<FnvHashSet<_>>();
        let expect = tets.into_iter().collect::<Vec<_>>();
        assert_eq!(result.len(), expect.len());
        let expect = expect
            .into_iter()
            .map(|[v0, v1, v2, v3]| even_sort_4([v(v0), v(v1), v(v2), v(v3)]))
            .collect::<FnvHashSet<_>>();
        assert_eq!(result, expect);
    }

    #[track_caller]
    fn assert_tets_ids<V, T, I: IntoIterator<Item = [VertexId; 4]>>(mesh: &TetMesh<V, T>, tets: I) {
        let result = mesh
            .tets
            .iter()
            .map(|(_, tet)| even_sort_4(tet.vertices()))
            .collect::<FnvHashSet<_>>();
        let expect = tets.into_iter().collect::<Vec<_>>();
        assert_eq!(result.len(), expect.len());
        let expect = expect
            .into_iter()
            .map(|[v0, v1, v2, v3]| even_sort_4([v0, v1, v2, v3]))
            .collect::<FnvHashSet<_>>();
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
        assert_eq!(
            mesh.vertex(TetMesh::<(), ()>::GHOST).map(|v| v.value()),
            None
        );

        assert_eq!(mesh.num_tets(), 5);
        assert_eq!(even_sort_4(mesh[t(0)].vertices()), [v(0), v(1), v(2), v(3)]);
        assert!(mesh[t(1)].vertices().contains(&TetMesh::<(), ()>::GHOST));
    }

    #[test]
    fn test_center() {
        let mesh = default_tets();
        assert_eq!(mesh.tet_volume_x6(t(0)), Some(1.0));
        assert_eq!(mesh.tet_centroid(t(0)), Some(Pt3::new(0.25, 0.25, 0.25)));
        assert_eq!(mesh.tet_volume_x6(t(1)), None);
        assert_eq!(mesh.tet_centroid(t(1)), None);
        assert_eq!(mesh.tet_volume_x6(t(2)), None);
        assert_eq!(mesh.tet_centroid(t(2)), None);
        assert_eq!(mesh.tet_volume_x6(t(3)), None);
        assert_eq!(mesh.tet_centroid(t(3)), None);
        assert_eq!(mesh.tet_volume_x6(t(4)), None);
        assert_eq!(mesh.tet_centroid(t(4)), None);
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
        assert_eq!(
            mesh.vertex(TetMesh::<(), ()>::GHOST).map(|v| v.value()),
            None
        );

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
        let walker = mesh.walker_from_tet_id(t(0));

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
        let walkers = mesh.walker_from_tet_id(t(1)).flip14_unchecked(&mut mesh, v(4));
        for (_, walker) in walkers.iter().enumerate() {
            assert_eq!(walker.fourth(&mesh), v(4));
        }

        mesh.assert_integrity(false);
    }

    #[test]
    fn test_flip23_unchecked() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let walkers = mesh.walker_from_tet_id(t(0)).flip14_unchecked(&mut mesh, v(4)); // ids 0, 5, 6, 7
        let walkers = walkers[0].flip23_unchecked(&mut mesh); // flips tri [2, 1, 3] away, opp verts 4 and GHOST
        assert_eq!(walkers[0].opp_edge(&mesh), walkers[1].opp_edge(&mesh));
        assert_eq!(walkers[1].opp_edge(&mesh), walkers[2].opp_edge(&mesh));
        assert!(walkers[0].opp_edge(&mesh).contains(&v(4)));
        assert!(walkers[0]
            .opp_edge(&mesh)
            .contains(&TetMesh::<(), ()>::GHOST));

        mesh.assert_integrity(false);
    }

    #[test]
    fn test_flip32_unchecked() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        let walkers = mesh.walker_from_tet_id(t(0)).flip14_unchecked(&mut mesh, v(4)); // ids 0, 5, 6, 7
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

        mesh.assert_integrity(false);
    }

    //#[test]
    //fn test_boundary_and_enclosure_all() {
    //    let mut mesh = default_tets();
    //    // Just care about combinatorial structure for this test
    //    mesh.vertices.insert(Vertex::new(
    //        TetWalker::new(TetId::invalid(), 0),
    //        Pt3::origin(),
    //        (),
    //    ));
    //    let (boundary, enclosed) = mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4))[0]
    //        .point_cavity(&mesh, |_| true);

    //    assert_eq!(boundary, vec![]);
    //    assert_eq!(
    //        enclosed.into_iter().collect::<FnvHashSet<_>>(),
    //        mesh.tets.keys().collect::<FnvHashSet<_>>()
    //    );
    //}

    //#[test]
    //fn test_boundary_and_enclosure_one() {
    //    let mut mesh = default_tets();
    //    // Just care about combinatorial structure for this test
    //    mesh.vertices.insert(Vertex::new(
    //        TetWalker::new(TetId::invalid(), 0),
    //        Pt3::origin(),
    //        (),
    //    ));

    //    let walker = mesh.walker_from_tet(t(0));
    //    let (boundary, enclosed) = walker.flip14_unchecked(&mut mesh, v(4))[0]
    //        .point_cavity(&mesh, |t| {
    //            mesh[t].vertices() == mesh[walker.id()].vertices()
    //        });

    //    assert_eq!(
    //        boundary
    //            .into_iter()
    //            .map(BoundaryTri)
    //            .collect::<FnvHashSet<_>>(),
    //        vec![
    //            BoundaryTri(walker),
    //            BoundaryTri(walker.to_twin_edge()),
    //            BoundaryTri(walker.to_opp_edge()),
    //            BoundaryTri(walker.to_perm(Permutation::_3210)),
    //        ]
    //        .into_iter()
    //        .collect::<FnvHashSet<_>>()
    //    );

    //    assert_eq!(enclosed, vec![walker.id()]);
    //}

    //#[test]
    //fn test_boundary_and_enclosure_some() {
    //    let mut mesh = default_tets();
    //    // Just care about combinatorial structure for this test
    //    mesh.vertices.insert(Vertex::new(
    //        TetWalker::new(TetId::invalid(), 0),
    //        Pt3::origin(),
    //        (),
    //    ));

    //    let walker = mesh.walker_from_tet(t(0));
    //    let (boundary, enclosed) = walker.flip14_unchecked(&mut mesh, v(4))[0]
    //        .point_cavity(&mesh, |t| {
    //            !mesh[t].vertices().contains(&TetMesh::<(), ()>::GHOST)
    //        });

    //    assert_eq!(
    //        boundary
    //            .into_iter()
    //            .map(|w| even_sort_3(w.tri(&mesh)))
    //            .collect::<FnvHashSet<_>>(),
    //        vec![
    //            [v(0), v(1), v(2)],
    //            [v(1), v(3), v(2)],
    //            [v(0), v(2), v(3)],
    //            [v(0), v(3), v(1)],
    //        ]
    //        .into_iter()
    //        .collect::<FnvHashSet<_>>()
    //    );

    //    assert_eq!(
    //        enclosed.into_iter().collect::<FnvHashSet<_>>(),
    //        vec![t(0), t(5), t(6), t(7)]
    //            .into_iter()
    //            .collect::<FnvHashSet<_>>()
    //    );
    //}

    //#[test]
    //fn test_flip_in_vertex_unchecked_1_4() {
    //    let mut mesh = default_tets();
    //    mesh.vertices.insert(Vertex::new(
    //        TetWalker::new(TetId::invalid(), 0),
    //        Pt3::origin(),
    //        (),
    //    ));
    //    let (mut boundary, enclosed) = mesh
    //        .walker_from_tet(t(1))
    //        .point_cavity(&mesh, |id| id == t(1));
    //    mesh.flip_in_vertex_unchecked(v(4), &mut boundary, &enclosed);

    //    for tet in 0..8 {
    //        for i in 0..12 {
    //            let walker = TetWalker::new(t(tet), i);
    //            let mut vertices = walker.tet(&mesh);

    //            if !vertices[..3].contains(&v(0)) {
    //                vertices.swap(0, 1);

    //                if vertices[3] == v(4) {
    //                    // Internal to external link
    //                    vertices[3] = v(0);
    //                } else if vertices[3] == v(0) {
    //                    // External to internal link
    //                    vertices[3] = v(4);
    //                } else {
    //                    // Internal link
    //                    vertices[3] = v((TetMesh::<(), ()>::GHOST.0 as u64 + 10
    //                        - vertices.iter().map(|v| v.0 as u64).sum::<u64>())
    //                        as IdType);
    //                }
    //                assert_eq!(
    //                    walker.to_adj(&mesh).tet(&mesh),
    //                    vertices,
    //                    "Flipping {:?}, expected {:?}, got {:?}",
    //                    walker.tet(&mesh),
    //                    vertices,
    //                    walker.to_adj(&mesh).tet(&mesh)
    //                );
    //            }
    //        }
    //    }

    //    mesh.assert_integrity(false);
    //}

    //#[test]
    //fn test_flip_in_vertex_unchecked_2_6() {
    //    let mut mesh = default_tets();
    //    for _ in 0..2 {
    //        mesh.vertices.insert(Vertex::new(
    //            TetWalker::new(TetId::invalid(), 0),
    //            Pt3::origin(),
    //            (),
    //        ));
    //    }

    //    mesh.walker_from_tet(t(0)).flip14_unchecked(&mut mesh, v(4));
    //    let (mut boundary, enclosed) =
    //        mesh.walker_from_tet(t(1))
    //            .point_cavity(&mesh, |t| {
    //                let vs = mesh[t].vertices();
    //                vs.contains(&v(1)) && vs.contains(&v(2)) && vs.contains(&v(3))
    //            });

    //    // 2-to-6 flip.
    //    mesh.flip_in_vertex_unchecked(v(5), &mut boundary, &enclosed);

    //    for tet in 0..12 {
    //        for i in 0..12 {
    //            let walker = TetWalker::new(t(tet), i);
    //            let mut vertices = walker.tet(&mesh);

    //            // Resultant tet or opposite resultant tet
    //            if vertices[3] == v(0) || vertices.contains(&v(5)) {
    //                vertices.swap(0, 1);

    //                if vertices[3] == v(5) {
    //                    // Internal to external link
    //                    vertices[3] = v(0);
    //                } else if vertices[3] == v(0) {
    //                    // External to internal link
    //                    vertices[3] = v(5);
    //                } else {
    //                    // Internal link
    //                    vertices[3] = if vertices[3] == v(4) {
    //                        TetMesh::<(), ()>::GHOST
    //                    } else if vertices[3] == TetMesh::<(), ()>::GHOST {
    //                        v(4)
    //                    } else {
    //                        v(6 - vertices
    //                            .iter()
    //                            .map(|v| v.0)
    //                            .filter(|n| *n == 1 || *n == 2 || *n == 3)
    //                            .sum::<u32>())
    //                    };
    //                }
    //                assert_eq!(
    //                    walker.to_adj(&mesh).tet(&mesh),
    //                    vertices,
    //                    "Flipping {:?}, expected {:?}, got {:?}",
    //                    walker.tet(&mesh),
    //                    vertices,
    //                    walker.to_adj(&mesh).tet(&mesh)
    //                );
    //            }
    //        }
    //    }

    //    mesh.assert_integrity(false);
    //}

    #[test]
    fn test_walkers_from_vertex() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        mesh.walker_from_tet_id(t(0)).flip14_unchecked(&mut mesh, v(4));

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
        mesh.walker_from_tet_id(t(0)).flip14_unchecked(&mut mesh, v(4));

        let tets = mesh.vertex_tets(v(4)).collect::<Vec<_>>();
        assert_eq!(tets.len(), 4);
        assert_eq!(
            tets.into_iter()
                .map(|tet| even_sort_4(mesh[tet].vertices()))
                .collect::<FnvHashSet<_>>(),
            vec![
                [v(1), v(2), v(4), v(3)],
                [v(0), v(2), v(3), v(4)],
                [v(0), v(1), v(4), v(3)],
                [v(0), v(1), v(2), v(4)],
            ]
            .into_iter()
            .collect::<FnvHashSet<_>>()
        );
    }

    #[test]
    fn test_vertex_targets() {
        let mut mesh = default_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::origin(),
            (),
        ));
        mesh.walker_from_tet_id(t(0)).flip14_unchecked(&mut mesh, v(4));

        let targets = mesh.vertex_targets(v(4)).collect::<Vec<_>>();
        assert_eq!(targets.len(), 4);
        assert_eq!(
            targets.into_iter().collect::<FnvHashSet<_>>(),
            vec![v(0), v(1), v(2), v(3)]
                .into_iter()
                .collect::<FnvHashSet<_>>()
        );
    }

    #[test]
    fn test_in_sphere_ghost() {
        let mut mesh = reverse_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::new(0.5, 0.3, 0.6),
            (),
        ));

        assert!(mesh.in_sphere_sim(TetMesh::<(), ()>::GHOST, v(1), v(2), v(3), v(4)));
        assert!(!mesh.in_sphere_sim(v(0), TetMesh::<(), ()>::GHOST, v(2), v(3), v(4)));
        assert!(!mesh.in_sphere_sim(v(0), v(1), TetMesh::<(), ()>::GHOST, v(3), v(4)));
        assert!(!mesh.in_sphere_sim(v(0), v(1), v(2), TetMesh::<(), ()>::GHOST, v(4)));
    }

    #[test]
    fn test_in_sphere_no_ghost() {
        let mut mesh = reverse_tets();
        mesh.vertices.insert(Vertex::new(
            TetWalker::new(TetId::invalid(), 0),
            Pt3::new(0.5, 0.3, 0.6),
            (),
        ));

        assert!(mesh.in_sphere_sim(v(0), v(1), v(3), v(2), v(4)));
    }

    #[test]
    fn test_delaunay_tets_single() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
            ],
            0.0,
            || (),
        );

        assert_eq!(mesh.num_vertices(), 4);
        assert_tets(
            &mesh,
            vec![
                [0, 1, 3, 2],
                [3, 1, 0, u32::MAX],
                [1, 3, 2, u32::MAX],
                [0, 2, 3, u32::MAX],
                [2, 0, 1, u32::MAX],
            ],
        );
    }

    #[test]
    fn test_delaunay_tets_multiple() {
        let mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(0.0, 0.0, 0.0), ()),
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(1.5, 1.5, 1.0), ()),
                (Pt3::new(0.5, 0.5, 0.5), ()),
            ],
            0.0,
            || (),
        );

        assert_eq!(mesh.num_vertices(), 6);
        assert_tets(
            &mesh,
            vec![
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
            ],
        );
    }

    #[test]
    fn test_progress_remove_edge_simple() {
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
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
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_edge([v(3), v(4)]).unwrap();
        assert_eq!(walker.edge(&mesh), [v(3), v(4)]);

        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        assert_tets(
            &mesh,
            vec![
                [0, 1, 2, 4],
                [0, 2, 1, 3],
                [0, 3, 1, u32::MAX],
                [1, 3, 2, u32::MAX],
                [2, 3, 0, u32::MAX],
                [0, 4, 2, u32::MAX],
                [2, 4, 1, u32::MAX],
                [1, 4, 0, u32::MAX],
            ],
        );
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_44_low_depth() {
        // Effectively a 4-to-4 flip, but fails due to the low depth.

        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 0.2), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_edge([v(4), v(5)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets(
            &mesh,
            vec![
                [0, 1, 4, 5],
                [1, 2, 4, 5],
                [2, 3, 4, 5],
                [3, 0, 4, 5],
                [1, 0, 4, u32::MAX],
                [2, 1, 4, u32::MAX],
                [3, 2, 4, u32::MAX],
                [0, 3, 4, u32::MAX],
                [0, 1, 5, u32::MAX],
                [1, 2, 5, u32::MAX],
                [2, 3, 5, u32::MAX],
                [3, 0, 5, u32::MAX],
            ],
        );
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_44() {
        // Effectively a 4-to-4 flip, but fails due to the low depth.

        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, 0.0, 0.0), ()),
                (Pt3::new(0.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 0.2), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_edge([v(4), v(5)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                2,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        assert_eq!(mesh.num_tets(), 12);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_concavity() {
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, 1.0), ()),
                (Pt3::new(0.0, 0.0, 0.1), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();
        let walker = mesh.walker_from_edge([v(3), v(4)]).unwrap();

        // Depth is high to try so hard and not get far at all.
        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                5,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_ghost() {
        // Successful removal of solid edge attached to ghost triangle
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(-0.8, 0.0, 0.2), ()),
                (Pt3::new(-0.8, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_edge([v(1), v(2)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        assert_tets(
            &mesh,
            vec![
                [0, 1, 3, 4],
                [3, 2, 0, 4],
                [3, 1, 0, u32::MAX],
                [0, 2, 3, u32::MAX],
                [4, 1, 3, u32::MAX],
                [3, 2, 4, u32::MAX],
                [0, 1, 4, u32::MAX],
                [4, 2, 0, u32::MAX],
            ],
        );
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_ghost_edge() {
        // Successful removal of ghost edge.
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(0.0, 0.0, -0.2), ()),
                (Pt3::new(0.0, 0.0, -2.0), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();

        let walker = mesh.walker_from_tri([v(0), v(1), v(2)]).unwrap();
        assert_eq!(walker.tri(&mesh), [v(0), v(1), v(2)]);
        walker.flip23_unchecked(&mut mesh); // Set up the ghost edge

        let walker = mesh
            .walker_from_edge([v(3), TetMesh::<(), ()>::GHOST])
            .unwrap();

        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        // Expect the Delaunay tetrahedralization
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_edge_ghost_concavity() {
        // An edge attached to a ghost triangle is concave.
        let mut mesh = default_tets();
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();
        let walker = mesh.walker_from_edge([v(0), v(1)]).unwrap();

        // Depth is high to try so hard and not get far at all.
        assert_eq!(
            walker.progress_by_removing_edge(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                5,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_simple() {
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
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
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_tri([v(0), v(1), v(2)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        assert_tets(
            &mesh,
            vec![
                [0, 1, 3, 4],
                [1, 2, 3, 4],
                [2, 0, 3, 4],
                [0, 3, 1, u32::MAX],
                [1, 3, 2, u32::MAX],
                [2, 3, 0, u32::MAX],
                [1, 4, 0, u32::MAX],
                [2, 4, 1, u32::MAX],
                [0, 4, 2, u32::MAX],
            ],
        );
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_complex_low_depth() {
        // Must remove an edge first, but depth is too low.

        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
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
        mesh.set_track_tri_map(true);
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();
        let walker = mesh.walker_from_tri([v(0), v(3), v(4)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_complex() {
        // Must remove an edge first.

        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
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
        mesh.set_track_tri_map(true);
        let walker = mesh.walker_from_tri([v(0), v(3), v(4)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                2,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        assert_tets(&mesh, vec![
            [0, 1, 2, 4],
            [0, 2, 1, 3],
            [0, 3, 1, u32::MAX],
            [1, 3, 2, u32::MAX],
            [2, 3, 0, u32::MAX],
            [0, 4, 2, u32::MAX],
            [2, 4, 1, u32::MAX],
            [1, 4, 0, u32::MAX],
        ]);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_concavity() {
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, -0.5, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(-1.1, 0.0, 2.0), ()),
                (Pt3::new(-1.1, 0.0, -2.0), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();
        let walker = mesh.walker_from_tri([v(0), v(1), v(2)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                5,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_ghost_concavity() {
        // An edge attached to a ghost triangle is concave.
        let mut mesh = default_tets();
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();
        let walker = mesh.walker_from_tri([v(0), v(1), v(2)]).unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                5,
                true,
                &mut vec![],
                &mut vec![]
            ),
            None
        );
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }

    #[test]
    fn test_progress_remove_tri_ghost_tri() {
        // Successful removal of ghost tri.
        let mut mesh = TetMesh::<(), ()>::delaunay_from_vertices(
            vec![
                (Pt3::new(1.0, 0.0, 0.0), ()),
                (Pt3::new(-1.0, 1.0, 0.0), ()),
                (Pt3::new(-1.0, -1.0, 0.0), ()),
                (Pt3::new(-0.6, 0.0, 0.2), ()),
                (Pt3::new(-0.6, 0.0, -0.2), ()),
            ],
            0.0,
            || (),
        );
        mesh.set_track_tri_map(true);
        let tets = mesh
            .tets()
            .map(|(_, tet)| tet.vertices())
            .collect::<Vec<_>>();

        mesh.walker_from_edge([v(1), v(2)]).unwrap().flip32_unchecked(&mut mesh);

        let walker = mesh
            .walker_from_tri([v(3), v(4), TetMesh::<(), ()>::GHOST])
            .unwrap();

        assert_eq!(
            walker.progress_by_removing_tri(
                &mut mesh,
                |_, _| true,
                |_, _| true,
                1,
                true,
                &mut vec![],
                &mut vec![]
            ),
            Some(0)
        );
        // Expect the Delaunay tetrahedralization
        assert_tets_ids(&mesh, tets);
        mesh.assert_integrity(true);
    }
}
