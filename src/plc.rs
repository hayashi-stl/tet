use crate::{Pt3, VertexId};
use crate::id_map::{self, IdMap};

crate::id! {
    /// A PLC edge id
    pub struct EdgeId
}

crate::id! {
    /// An id for a specific edge of a PLC face
    pub struct FaceEdgeId
}

crate::id! {
    /// An id for a ring of edges around a PLC face
    pub struct FaceRingId
}

crate::id! {
    /// A PLC face id
    pub struct FaceId
}

/// A PLC vertex.
#[derive(Clone, Debug)]
pub struct Vertex<V> {
    /// An edge that it's the source of
    edge_out: EdgeId,
    /// An edge that it's the target of
    edge_in: EdgeId,
    position: Pt3,
    value: V,
}

impl<V> Vertex<V> {
    fn new(position: Pt3, value: V) -> Self {
        Self {
            edge_out: EdgeId::invalid(),
            edge_in: EdgeId::invalid(),
            position,
            value,
        }
    }

    fn position(&self) -> Pt3 {
        self.position
    }

    fn value(&self) -> &V {
        &self.value
    }
}

/// A PLC edge.
/// The vertices are the same if this edge is an isolated vertex on a face.
#[derive(Clone, Debug)]
pub struct Edge<E> {
    vertices: [VertexId; 2],
    next_edge_out: EdgeId,
    next_edge_in: EdgeId,
    prev_edge_out: EdgeId,
    prev_edge_in: EdgeId,
    /// A face edge that it's part of
    face_edge: FaceEdgeId,
    value: E,
}

impl<E> Edge<E> {
    fn new(vertices: [VertexId; 2], value: E) -> Self {
        Self {
            vertices,
            next_edge_out: EdgeId::invalid(),
            next_edge_in: EdgeId::invalid(),
            prev_edge_out: EdgeId::invalid(),
            prev_edge_in: EdgeId::invalid(),
            face_edge: FaceEdgeId::invalid(),
            value,
        }
    }

    fn vertices(&self) -> [VertexId; 2] {
        self.vertices
    }

    fn value(&self) -> &E {
        &self.value
    }
}

/// An edge of a PLC face.
#[derive(Clone, Debug)]
pub struct FaceEdge {
    edge: EdgeId,
    ring: FaceRingId,
    /// Next face edge along face
    next: FaceEdgeId,
    /// Previous face edge along face
    prev: FaceEdgeId,
    /// Next face edge on this edge
    next_on_edge: FaceEdgeId,
    /// Prev face edge on this edge
    prev_on_edge: FaceEdgeId,
}

impl FaceEdge {
    fn new(edge: EdgeId, ring: FaceRingId) -> Self {
        Self {
            edge,
            ring,
            next: FaceEdgeId::invalid(),
            prev: FaceEdgeId::invalid(),
            next_on_edge: FaceEdgeId::invalid(),
            prev_on_edge: FaceEdgeId::invalid(),
        }
    }
}

/// A ring of edges around a PLC face.
#[derive(Clone, Debug)]
pub struct FaceRing {
    /// A face edge that belongs to this ring
    face_edge: FaceEdgeId,
    face: FaceId,
    next: FaceRingId,
    prev: FaceRingId,
}

impl FaceRing {
    fn new(face: FaceId) -> Self {
        Self {
            face_edge: FaceEdgeId::invalid(),
            face,
            next: FaceRingId::invalid(),
            prev: FaceRingId::invalid(),
        }
    }
}

/// A PLC face.
#[derive(Clone, Debug)]
pub struct Face<F> {
    /// A face ring that has this as its face
    ring: FaceRingId,
    value: F,
}

impl<F> Face<F> {
    fn new(value: F) -> Self {
        Self {
            ring: FaceRingId::invalid(),
            value,
        }
    }

    fn value(&self) -> &F {
        &self.value
    }
}

/// A piecewise linear complex. Contains vertices, edges and faces.
/// The faces do not have to be triangles. In fact,
/// they can contain holes, slits, and internal vertices.
#[derive(Clone, Debug)]
pub struct Plc<V, E, F> {
    vertices: IdMap<VertexId, Vertex<V>>,
    edges: IdMap<EdgeId, Edge<E>>,
    face_edges: IdMap<FaceEdgeId, FaceEdge>,
    rings: IdMap<FaceRingId, FaceRing>,
    faces: IdMap<FaceId, Face<F>>,
    default_vertex: fn() -> V,
    default_edge: fn() -> E,
    default_face: fn() -> F,
}

impl<V, E, F> Plc<V, E, F> {
    /// Creates an empty PLC.
    pub fn new(default_vertex: fn() -> V, default_edge: fn() -> E, default_face: fn() -> F) -> Self {
        Self {
            vertices: IdMap::default(),
            edges: IdMap::default(),
            face_edges: IdMap::default(),
            rings: IdMap::default(),
            faces: IdMap::default(),
            default_vertex,
            default_edge,
            default_face,
        }
    }

    /// Gets the number of vertices
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets the number of edges
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Gets the number of faces
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    /// Gets a vertex, if it exists.
    pub fn vertex(&self, vertex: VertexId) -> Option<&Vertex<V>> {
        self.vertices.get(vertex)
    }

    /// Gets a edge, if it exists.
    pub fn edge(&self, edge: EdgeId) -> Option<&Edge<E>> {
        self.edges.get(edge)
    }

    /// Gets a face, if it exists.
    pub fn face(&self, face: FaceId) -> Option<&Face<F>> {
        self.faces.get(face)
    }

    /// Iterates over vertex ids and vertices.
    pub fn vertices(&self) -> Vertices<V> {
        self.vertices.iter()
    }

    /// Iterates over edge ids and edges.
    pub fn edges(&self) -> Edges<E> {
        self.edges.iter()
    }

    /// Iterates over face ids and faces.
    pub fn faces(&self) -> Faces<F> {
        self.faces.iter()
    }
}

/// Iterates over vertex ids and vertices.
pub type Vertices<'a, V> = id_map::Iter<'a, VertexId, Vertex<V>>;

/// Iterates over edge ids and edges.
pub type Edges<'a, E> = id_map::Iter<'a, EdgeId, Edge<E>>;

/// Iterates over face ids and faces.
pub type Faces<'a, F> = id_map::Iter<'a, FaceId, Face<F>>;