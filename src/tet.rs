use crate::{VertexId, TetId};

/// A vertex stores a tet that it's part of.
struct Vertex<V> {
    tet: TetId,
    value: V,
}

/// A tet stores the vertices that it's part of (in the correct orientation)
/// and the opposite tets of each vertex.
///
/// Walkers store the index of the first edge of the opposite triangle.
struct Tet<T> {
    vertices: [VertexId; 4],
    opp_tets: [TetWalker; 4],
    value: T,
}

/// A walker over the tetrahedrons of a tet mesh.
/// The edge index is a number in 0..12 indexing into
/// the edges of a tet.
///
/// The enumeration of the edges is borrowed from TetGen.
/// One edge of each triangle is enumerated,
/// then the next edge of each triangle respectively,
/// then the last edge of each triangle respectively.
struct TetWalker {
    tet: TetId,
    edge_index: u32,
}