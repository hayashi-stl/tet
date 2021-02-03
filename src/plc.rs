use crate::{Pt3, VertexId};
use crate::id_map::{self, IdMap};

use fnv::FnvHashMap;
use std::{collections::hash_map::Entry, iter};

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
    /// The face ring that it basically is if it's an isolated vertex of a face.
    ring: FaceRingId,
    position: Pt3,
    value: V,
}

impl<V> Vertex<V> {
    fn new(position: Pt3, value: V) -> Self {
        Self {
            edge_out: EdgeId::invalid(),
            edge_in: EdgeId::invalid(),
            ring: FaceRingId::invalid(),
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Element {
    FaceEdge(FaceEdgeId),
    Vertex(VertexId),
}

/// A ring of edges around a PLC face.
#[derive(Clone, Debug)]
pub struct FaceRing {
    /// A face edge or vertex that belongs to this ring.
    /// This represents an isolated vertex on a face if it's a vertex.
    element: Element,
    face: FaceId,
    next: FaceRingId,
    prev: FaceRingId,
}

impl FaceRing {
    fn new(face: FaceId) -> Self {
        Self {
            element: Element::FaceEdge(FaceEdgeId::invalid()),
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
    edge_map: FnvHashMap<[VertexId; 2], EdgeId>,
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
            edge_map: FnvHashMap::default(),
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

    /// Adds a vertex to the PLC and returns its vertex id.
    pub fn add_vertex(&mut self, position: Pt3, value: V) -> VertexId {
        self.vertices.insert(Vertex::new(position, value))
    }

    /// Adds an edge to the PLC and returns its edge id, assuming the edge doesn't already exist.
    /// Does not add it to the hash map; that should be done manually.
    fn add_new_edge(vertex_list: &mut IdMap<VertexId, Vertex<V>>, edge_list: &mut IdMap<EdgeId, Edge<E>>,
        vertices: [VertexId; 2], value: E) -> EdgeId
    {
        let id = edge_list.insert(Edge::new(vertices, value));

        // Fix source links
        let out_id = vertex_list[vertices[0]].edge_out;
        if out_id == EdgeId::invalid() {
            vertex_list[vertices[0]].edge_out = id;
            edge_list[id].next_edge_out = id;
            edge_list[id].prev_edge_out = id;
        } else {
            let prev = edge_list[out_id].prev_edge_out;
            edge_list[prev].next_edge_out = id;
            edge_list[out_id].prev_edge_out = id;
            edge_list[id].next_edge_out = out_id;
            edge_list[id].prev_edge_out = out_id;
        }

        // Fix target links
        let in_id = vertex_list[vertices[1]].edge_in;
        if in_id == EdgeId::invalid() {
            vertex_list[vertices[1]].edge_in = id;
            edge_list[id].next_edge_in = id;
            edge_list[id].prev_edge_in = id;
        } else {
            let prev = edge_list[in_id].prev_edge_in;
            edge_list[prev].next_edge_in = id;
            edge_list[in_id].prev_edge_in = id;
            edge_list[id].next_edge_in = in_id;
            edge_list[id].prev_edge_in = in_id;
        }

        id
    }

    /// Adds an edge to the PLC and returns its edge id.
    /// Simply sets the value if the edge already exists.
    pub fn add_edge(&mut self, vertices: [VertexId; 2], value: E) -> EdgeId {
        match self.edge_map.entry(vertices) {
            Entry::Vacant(entry) => {
                let id = Self::add_new_edge(&mut self.vertices, &mut self.edges, vertices, value);
                entry.insert(id);
                id
            }

            Entry::Occupied(entry) => {
                let id = *entry.get();
                self.edges[id].value = value;
                id
            }
        }
    }

    /// Adds an edge to the PLC and returns its edge id.
    /// Has no side effects if the edge already exists.
    pub fn add_edge_if_absent(&mut self, vertices: [VertexId; 2], value: E) -> EdgeId {
        match self.edge_map.entry(vertices) {
            Entry::Vacant(entry) => {
                let id = Self::add_new_edge(&mut self.vertices, &mut self.edges, vertices, value);
                entry.insert(id);
                id
            }

            Entry::Occupied(entry) => {
                *entry.get()
            }
        }
    }
    
    /// Adds a face to the PLC and returns its id.
    ///
    /// Warning: This can cause 2 faces with the same vertices in the same order to exist.
    pub fn add_face<RI: IntoIterator<Item = VI>, VI: IntoIterator<Item = VertexId>>(&mut self, vertices: RI, value: F) -> FaceId {
        let id = self.faces.insert(Face::new(value));
        
        // Rings
        let mut first_ring = FaceRingId::invalid();
        let mut prev_ring = FaceRingId::invalid();

        for ring_iter in vertices.into_iter() {
            let ring = self.rings.insert(FaceRing::new(id));
            if first_ring == FaceRingId::invalid() {
                first_ring = ring;
                self.faces[id].ring = first_ring;
            } else {
                self.rings[prev_ring].next = ring;
                self.rings[ring].prev = prev_ring;
            }

            // Face edges
            let mut ring_iter = ring_iter.into_iter();
            let mut first_vertex = ring_iter.next().expect("Face ring must have at least 1 vertex.");
            let mut prev_vertex = first_vertex;

            if let Some(second_vertex) = ring_iter.next() {
                let mut first_fe = FaceEdgeId::invalid();
                let mut prev_fe = FaceEdgeId::invalid();

                // Multiple vertices
                for vertex in iter::once(second_vertex).chain(ring_iter).chain(iter::once(first_vertex)) {
                    let edge = self.add_edge_if_absent([prev_vertex, vertex], (self.default_edge)());
                    let fe = self.face_edges.insert(FaceEdge::new(edge, ring));

                    if first_fe == FaceEdgeId::invalid() {
                        first_fe = fe;
                        self.rings[ring].element = Element::FaceEdge(first_fe);
                    } else {
                        self.face_edges[prev_fe].next = fe;
                        self.face_edges[fe].prev = prev_fe;
                    }

                    // TODO: Edge links

                    prev_vertex = vertex;
                    prev_fe = fe;
                }

                // Complete cycle
                self.face_edges[prev_fe].next = first_fe;
                self.face_edges[first_fe].prev = prev_fe;
            } else {
                // Isolated vertex
                self.rings[ring].element = Element::Vertex(first_vertex);
            }

            prev_ring = ring;
        }

        // Complete cycle
        self.rings[prev_ring].next = first_ring;
        self.rings[first_ring].prev = prev_ring;

        id
    }
}

/// Iterates over vertex ids and vertices.
pub type Vertices<'a, V> = id_map::Iter<'a, VertexId, Vertex<V>>;

/// Iterates over edge ids and edges.
pub type Edges<'a, E> = id_map::Iter<'a, EdgeId, Edge<E>>;

/// Iterates over face ids and faces.
pub type Faces<'a, F> = id_map::Iter<'a, FaceId, Face<F>>;

///// Iterates over vertex ids generated from extending the vertex list.
//pub type ExtendVertices<'a, V, I> = id_map::ExtendValues<'a, VertexId, Vertex<V>, I>;
//
///// Iterates over edge ids generated from extending the edge list.
//pub type ExtendEdges<'a, E, I> = id_map::ExtendValues<'a, EdgeId, Edge<E>, I>;
//
///// Iterates over face ids generated from extending the face list.
//pub type ExtendFaces<'a, F, I> = id_map::ExtendValues<'a, FaceId, Face<F>, I>;