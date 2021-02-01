use std::ops::Index;

use crate::id_map::{IdMap, IdType};
use crate::{Pt3, TetId, VertexId};
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
    edge_index: u32,
}

impl TetWalker {
    fn new(tet: TetId, edge_index: u32) -> Self {
        Self { tet, edge_index }
    }

    /// Gets the first vertex of the tet walker. This is the current vertex.
    pub fn first<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2][self.edge_index as usize]]
    }

    /// Gets the second vertex of the tet walker. This is the current edge's target.
    pub fn second<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[2, 3, 0, 1, 1, 0, 3, 2, 3, 2, 1, 0][self.edge_index as usize]]
    }

    /// Gets the third vertex of the tet walker.
    pub fn third<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1][self.edge_index as usize]]
    }

    /// Gets the fourth vertex of the tet walker. This is the vertex opposite the current triangle.
    pub fn fourth<V, T>(self, mesh: &Tets<V, T>) -> VertexId {
        mesh.tets[self.tet].vertices[[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3][self.edge_index as usize]]
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
            self.edge_index = (self.edge_index + 4) % 12;
            self
        }
    }

    crate::alias! {
        to_pfe,
        /// Advance to the previous edge in the current triangle.
        pub fn to_prev_tri_edge((mut) self: Self) -> Self {
            self.edge_index = (self.edge_index + 8) % 12;
            self
        }
    }

    crate::alias! {
        to_nve,
        /// Advance to the next edge from the current vertex. This is a counterclockwise rotation.
        pub fn to_next_vertex_edge((mut) self: Self) -> Self {
            self.edge_index = [10, 11, 8, 9, 1, 0, 3, 2, 7, 6, 5, 4][self.edge_index as usize];
            self
        }
    }

    crate::alias! {
        to_pve,
        /// Advance to the previous edge from the current vertex. This is a clockwise rotation.
        pub fn to_prev_vertex_edge((mut) self: Self) -> Self {
            self.edge_index = [5, 4, 7, 6, 11, 10, 9, 8, 2, 3, 0, 1][self.edge_index as usize];
            self
        }
    }

    /// Flip the current edge, keeping the current tet the same.
    pub fn to_twin_edge(mut self) -> Self {
        self.edge_index = [1, 0, 3, 2, 7, 6, 5, 4, 10, 11, 8, 9][self.edge_index as usize];
        self
    }

    /// Sets the current edge to the opposite edge of the current tet.
    pub fn to_opp_edge(mut self) -> Self {
        self.edge_index = [2, 3, 0, 1, 5, 4, 7, 6, 11, 10, 9, 8][self.edge_index as usize];
        self
    }

    crate::alias! {
        to_adj,
        /// Flip the current triangle, moving to the adjacent tet on that triangle.
        /// This flips the current edge.
        pub fn to_twin_tri<V, T>((mut) self: Self, mesh: &Tets<V, T>) -> Self {
            let div = self.edge_index / 4;
            self = mesh[self.id()].opp_tets[self.edge_index as usize % 4];
            match div {
                0 => self,
                1 => self.to_prev_tri_edge(),
                2 => self.to_next_tri_edge(),
                _ => unreachable!(),
            }
        }
    }

    /// Sets the orientation such that the vertices returned by `self.vertices()`
    /// are in the same order as those returned by `mesh[self.id()].vertices()`.
    pub fn to_canon_tet(mut self) -> Self {
        self.edge_index = 3;
        self
    }

    /// Performs a 1-to-4 flip on the mesh without checking whether
    /// such a flip produces negative-volume tets. Use at your own risk.
    /// Moves the walker to the new vertex while keeping the opposite triangle the same.
    /// Also returns the 4 tet ids from the flip, with the first one being reused.
    pub fn flip14_unchecked<V: Clone, T: Clone>(self, mesh: &mut Tets<V, T>, vertex: VertexId) -> (Self, [TetId; 4]) {
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
            mesh.tets[walker.id()].opp_tets[walker.edge_index as usize % 4].tet = ids[i];
        }

        (self, ids)
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
            .expect(&format!("Vertex {} does not exist", index))
    }
}

impl<V, T> Index<TetId> for Tets<V, T> {
    type Output = Tet<T>;

    fn index(&self, index: TetId) -> &Self::Output {
        self.tet(index)
            .expect(&format!("Tet {} does not exist", index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        mesh.vertices.insert(Vertex::new(TetWalker::new(TetId::invalid(), 0), Pt3::origin(), ()));
        mesh.walker_from_tet(t(1)).flip14_unchecked(&mut mesh, v(4));

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
    }
}
