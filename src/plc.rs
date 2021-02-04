mod intersect;

use crate::{Pt1, Pt3, Vec1, Vec3, VertexId, util::{CircularListIter, MapWith}};
use crate::id_map::{self, IdMap};

use float_ord::FloatOrd;
use fnv::{FnvHashMap, FnvHashSet};
use iter::Map;
use nalgebra::Unit;
use ncollide3d::partitioning::{BVH, BVT};
use std::{cell::Cell, collections::hash_map::Entry, fmt::Debug, iter, ops::Index};
use std::hash::Hash;
use std::fmt::{Display, Formatter, self};
use std::error::Error;

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

    /// Gets the position of this vertex.
    pub fn position(&self) -> Pt3 {
        self.position
    }

    /// Gets the value of this vertex.
    pub fn value(&self) -> &V {
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

    /// Gets the vertices of this edge.
    pub fn vertices(&self) -> [VertexId; 2] {
        self.vertices
    }

    /// Gets the value of this edge.
    pub fn value(&self) -> &E {
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
enum RingElement {
    FaceEdge(FaceEdgeId),
    Vertex(VertexId),
}

/// A ring of edges around a PLC face.
#[derive(Clone, Debug)]
pub struct FaceRing {
    /// A face edge or vertex that belongs to this ring.
    /// This represents an isolated vertex on a face if it's a vertex.
    element: RingElement,
    face: FaceId,
    next: FaceRingId,
    prev: FaceRingId,
}

impl FaceRing {
    fn new(face: FaceId) -> Self {
        Self {
            element: RingElement::FaceEdge(FaceEdgeId::invalid()),
            face,
            next: FaceRingId::invalid(),
            prev: FaceRingId::invalid(),
        }
    }

    /// Iterate over the vertices of this face ring.
    pub fn vertices<'a, V, E, F>(&self, plc: &'a Plc<V, E, F>) -> FaceRingVertices<'a, V, E, F> {
        match self.element {
            RingElement::FaceEdge(fe) =>
                FaceRingVertices::Multiple(
                    CircularListIter::new(plc, fe, |plc, fe| &plc.face_edges[fe], |plc, fe| plc.face_edges[fe].next)
                ),

            RingElement::Vertex(v) => FaceRingVertices::Single(Some((v, &plc.vertices[v])))
        }
    }

    /// Iterate over the edges of this face ring.
    pub fn edges<'a, V, E, F>(&self, plc: &'a Plc<V, E, F>) -> FaceRingEdges<'a, V, E, F> {
        match self.element {
            RingElement::FaceEdge(fe) =>
                FaceRingEdges::Multiple(
                    CircularListIter::new(plc, fe, |plc, fe| &plc.face_edges[fe], |plc, fe| plc.face_edges[fe].next)
                ),

            RingElement::Vertex(v) => FaceRingEdges::Single
        }
    }
}

/// A PLC face.
#[derive(Clone, Debug)]
pub struct Face<F> {
    /// A face ring that has this as its face
    ring: FaceRingId,
    /// The plane the face lies on in (point, normal) form.
    plane: (Pt3, Unit<Vec3>),
    value: F,
}

impl<F> Face<F> {
    fn new(value: F) -> Self {
        Self {
            ring: FaceRingId::invalid(),
            plane: (Pt1::new(f64::NAN).xxx(), Unit::new_unchecked(Vec3::new(1.0, 0.0, 0.0))),
            value,
        }
    }

    /// Gets the value of this face.
    pub fn value(&self) -> &F {
        &self.value
    }

    /// Iterates over the rings of this face.
    pub fn rings<'a, V, E>(&self, plc: &'a Plc<V, E, F>) -> FaceRings<'a, V, E, F> {
        CircularListIter::new(plc, self.ring, |plc, id| &plc.rings[id], |plc, id| plc.rings[id].next)
    }

    /// Iterates over the vertices of this face, one ring at a time.
    pub fn vertices<'a, V, E>(&self, plc: &'a Plc<V, E, F>) -> FaceVerticesByRing<'a, V, E, F> {
        FaceVerticesByRing(self.rings(plc))
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

    /// Iterates over vertex ids
    pub fn vertex_ids(&self) -> VertexIds<V> {
        self.vertices.keys()
    }

    /// Iterates over edge ids
    pub fn edge_ids(&self) -> EdgeIds<E> {
        self.edges.keys()
    }

    /// Iterates over face ids
    pub fn face_ids(&self) -> FaceIds<F> {
        self.faces.keys()
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
            edge_list[id].prev_edge_out = prev;
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
            edge_list[id].prev_edge_in = prev;
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
            let first_vertex = ring_iter.next().expect("Face ring must have at least 1 vertex.");
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
                        self.rings[ring].element = RingElement::FaceEdge(first_fe);
                    } else {
                        self.face_edges[prev_fe].next = fe;
                        self.face_edges[fe].prev = prev_fe;
                    }

                    let next_on_edge = self.edges[edge].face_edge;
                    if next_on_edge == FaceEdgeId::invalid() {
                        self.edges[edge].face_edge = fe;
                        self.face_edges[fe].next_on_edge = fe;
                        self.face_edges[fe].prev_on_edge = fe;
                    } else {
                        let prev_on_edge = self.face_edges[next_on_edge].prev_on_edge;
                        self.face_edges[next_on_edge].prev_on_edge = fe;
                        self.face_edges[prev_on_edge].next_on_edge = fe;
                        self.face_edges[fe].next_on_edge = next_on_edge;
                        self.face_edges[fe].prev_on_edge = prev_on_edge;
                    }

                    prev_vertex = vertex;
                    prev_fe = fe;
                }

                // Complete cycle
                self.face_edges[prev_fe].next = first_fe;
                self.face_edges[first_fe].prev = prev_fe;
            } else {
                // Isolated vertex
                self.rings[ring].element = RingElement::Vertex(first_vertex);
                self.vertices[first_vertex].ring = ring;
            }

            prev_ring = ring;
        }

        // Complete cycle
        self.rings[prev_ring].next = first_ring;
        self.rings[first_ring].prev = prev_ring;

        self.calc_face_plane(id);
        id
    }

    fn calc_face_plane(&mut self, face: FaceId) {
        let face_ = &self.faces[face];

        let points = face_.rings(self).flat_map(|(_, ring)| ring.vertices(self))
            .map(|(_, vertex)| vertex.position()).collect::<Vec<_>>();

        let center: Pt3 = (points.iter().map(|p| p.coords).sum::<Vec3>() / points.len() as f64).into();
        let normal = face_.rings(self).flat_map(|(_, ring)| ring.edges(self))
            .map(|(_, edge)| {
                let vs = edge.vertices();
                (self[vs[0]].position() - center).cross(&(self[vs[1]].position() - center))
            }).sum::<Vec3>();

        self.faces[face].plane = (center, Unit::new_normalize(normal))
    }

    /// Check that this is a valid PLC.
    /// That means all faces have a ring with at least 3 vertices,
    /// the vertices of each face are coplanar,
    /// and no two elements intersect each other without sharing a vertex.
    pub fn validate(&self, distance_tolerance: f64) -> Result<(), ValidateError<V, E, F>> {
        let mut degenerate = vec![];
        for (id, face) in self.faces() {
            if !face.rings(self).any(|(_, ring)| ring.vertices(self).count() >= 3) {
                degenerate.push(id);
            }
        }
        
        if !degenerate.is_empty() {
            Err(ValidateError::new(self, ValidateReason::DegenerateFace(degenerate), distance_tolerance))?;
        }

        let mut twisted = vec![];
        for (id, face) in self.faces() {
            let points = face.rings(self).flat_map(|(_, ring)| ring.vertices(self))
                .map(|(id, vertex)| (id, vertex.position())).collect::<Vec<_>>();

            let mut fails = vec![];
            for (vertex, point) in points {
                // Half the tolerance because the plane is in the middle.
                // Intersection tests should still work near the edge of the coplanarity tolerance.
                let dist = (point - face.plane.0).dot(&face.plane.1);
                if dist >= distance_tolerance / 2.0 {
                    fails.push((vertex, dist));
                    break;
                }
            }
            if !fails.is_empty() {
                twisted.push((id, fails));
            }
        }

        if !twisted.is_empty() {
            Err(ValidateError::new(self, ValidateReason::TwistedFace(twisted), distance_tolerance))?;
        }

        // Intersection time.
        let bvt = BVT::new_balanced(
            self.vertex_ids().map(|v| intersect::vertex_bound(self, v, distance_tolerance)).chain(
                self.edge_ids().map(|e| intersect::edge_bound(self, e, distance_tolerance))
            ).chain(self.face_ids().map(|f| intersect::face_bound(self, f, distance_tolerance))).collect()
        );

        let mut visitor = intersect::IntersectionTest {
            plc: self,
            tolerance: distance_tolerance,
            intersections: vec![],
        };

        bvt.visit_bvtt(&bvt, &mut visitor);
        if !visitor.intersections.is_empty() {
            Err(ValidateError::new(self, ValidateReason::Intersects(visitor.intersections), distance_tolerance))?;
        }

        Ok(())
    }

    #[cfg(feature = "obj")]
    pub fn from_obj(data: obj::ObjData, default_vertex: fn() -> V, default_edge: fn() -> E, default_face: fn() -> F) -> Self {
        let mut plc = Self::new(default_vertex, default_edge, default_face);

        for [x, y, z] in data.position {
            plc.add_vertex(Pt3::new(x as f64, y as f64, z as f64), default_vertex());
        }
        
        for object in data.objects {
            for group in object.groups {
                for obj::SimplePolygon(poly) in group.polys {
                    if poly.len() == 2 {
                        // Edge
                        plc.add_edge([VertexId(poly[0].0 as u32), VertexId(poly[1].0 as u32)], default_edge());
                    } else {
                        plc.add_face(
                            iter::once(poly.into_iter().map(|vertex| VertexId(vertex.0 as u32))),
                            default_face(),
                        );
                    }
                }
            }
        }

        plc
    }

    #[cfg(feature = "obj")]
    pub fn load_obj<P: AsRef<std::path::Path>>(path: P, default_vertex: fn() -> V, default_edge: fn() -> E, default_face: fn() -> F)
    -> Result<Self, obj::ObjError> {
        obj::Obj::load(path).map(|obj| Self::from_obj(obj.data, default_vertex, default_edge, default_face))
    }
}

/// An error returned from the validate function.
#[derive(Clone)]
pub struct ValidateError<'a, V, E, F> {
    /// Needed to display error
    plc: &'a Plc<V, E, F>,
    reason: ValidateReason,
    tolerance: f64,
}

impl<'a, V, E, F> ValidateError<'a, V, E, F> {
    fn new(plc: &'a Plc<V, E, F>, reason: ValidateReason, tolerance: f64) -> Self {
        Self { plc, reason, tolerance }
    }

    /// The reason why validation failed
    pub fn reason(&self) -> &ValidateReason {
        &self.reason
    }

    pub fn fmt_vertex(&self, vertex: VertexId, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "vertex {}", vertex)
    }

    pub fn fmt_edge(&self, edge: EdgeId, f: &mut Formatter<'_>) -> fmt::Result {
        let vs = self.plc[edge].vertices();
        write!(f, "edge {} ({}-{})", edge, vs[0], vs[1])
    }

    pub fn fmt_face(&self, face: FaceId, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "face {} [", face)?;
        let mut passed_first = false;

        for mut ring in self.plc[face].vertices(self.plc) {
            if passed_first {
                write!(f, ", ")?;
            } else {
                passed_first = true;
            }

            let first = ring.next().unwrap().0;
            write!(f, "{}", first)?;
            for vertex in ring.into_iter().map(|(v, _)| v).chain(iter::once(first)) {
                write!(f, "-{}", vertex)?;
            }
        }
        write!(f, "]")
    }

    pub fn fmt_element(&self, element: Element, f: &mut Formatter<'_>) -> fmt::Result {
        match element {
            Element::Vertex(vertex) => self.fmt_vertex(vertex, f),
            Element::Edge(edge) => self.fmt_edge(edge, f),
            Element::Face(face) => self.fmt_face(face, f),
        }
    }
}

impl<'a, V, E, F> Display for ValidateError<'a, V, E, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match &self.reason {
            ValidateReason::DegenerateFace(faces) => {
                writeln!(f, "The following faces have no ring with at least 3 vertices:")?;
                for face in faces {
                    self.fmt_face(*face, f)?;
                    writeln!(f)?
                }
            },

            ValidateReason::TwistedFace(faces) => {
                writeln!(f, "The following faces are twisted (not coplanar) (tolerance {}):", self.tolerance / 2.0)?;
                for (face, fails) in faces {
                    self.fmt_face(*face, f)?;
                    write!(f, " (failing vertices (distances from plane): {} ({})", fails[0].0, fails[0].1)?;
                    for (vertex, dist) in fails.into_iter().skip(1) {
                        write!(f, ", {} ({})", vertex, dist)?;
                    }
                    writeln!(f, ")")?;
                }
            }

            ValidateReason::Intersects(intersections) => {
                writeln!(f, "The following elements intersect or are too close together (tolerance {}):", self.tolerance)?;
                for (e0, e1, dist) in intersections {
                    self.fmt_element(*e0, f)?;
                    write!(f, " and ")?;
                    self.fmt_element(*e1, f)?;
                    write!(f, " (distance {})", dist)?;
                    writeln!(f)?
                }
            }
        };

        Ok(())
    }
}

impl<'a, V, E, F> Debug for ValidateError<'a, V, E, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Display::fmt(&self, f)
    }
}

impl<'a, V, E, F> Error for ValidateError<'a, V, E, F> {}

/// The error that actually happened.
#[derive(Clone, Debug, PartialEq)]
pub enum ValidateReason {
    /// Some faces have less than 3 vertices.
    DegenerateFace(Vec<FaceId>),
    /// Some faces are not coplanar enough.
    TwistedFace(Vec<(FaceId, Vec<(VertexId, f64)>)>),
    /// Some elements are too close to each other.
    Intersects(Vec<(Element, Element, f64)>),
}

/// A vertex, edge, or face.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Element {
    Vertex(VertexId),
    Edge(EdgeId),
    Face(FaceId),
}

/// Iterator over the rings of a face
pub type FaceRings<'a, V, E, F> = CircularListIter<&'a Plc<V, E, F>, FaceRingId, fn(&'a Plc<V, E, F>, FaceRingId) -> &'a FaceRing,
    fn(&'a Plc<V, E, F>, FaceRingId) -> FaceRingId>;

/// Iterator over the vertices of a face, one ring at a time
pub struct FaceVerticesByRing<'a, V, E, F>(FaceRings<'a, V, E, F>);

impl<'a, V, E, F> Iterator for FaceVerticesByRing<'a, V, E, F> {
    type Item = FaceRingVertices<'a, V, E, F>;

    fn next(&mut self) -> Option<Self::Item> {
        let plc = self.0.common;
        self.0.next().map(|(_, ring)| ring.vertices(plc))
    }
}

/// Iterator over the vertices of a face ring
pub enum FaceRingVertices<'a, V, E, F> {
    Single(Option<(VertexId, &'a Vertex<V>)>),
    Multiple(
        CircularListIter<&'a Plc<V, E, F>, FaceEdgeId, fn(&'a Plc<V, E, F>, FaceEdgeId) -> &'a FaceEdge, fn(&'a Plc<V, E, F>, FaceEdgeId) -> FaceEdgeId>,
    )
}

impl<'a, V, E, F> Iterator for FaceRingVertices<'a, V, E, F> {
    type Item = (VertexId, &'a Vertex<V>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Single(vertex) => vertex.take(),

            Self::Multiple(iter) => {
                let plc = iter.common;
                iter.next().map(|(_, fe)| {
                    let vertex = plc.edges[fe.edge].vertices()[0];
                    (vertex, &plc.vertices[vertex])
                })
            }
        }
    }
}

/// Iterator over the edges of a face ring
pub enum FaceRingEdges<'a, V, E, F> {
    Single,
    Multiple(
        CircularListIter<&'a Plc<V, E, F>, FaceEdgeId, fn(&'a Plc<V, E, F>, FaceEdgeId) -> &'a FaceEdge, fn(&'a Plc<V, E, F>, FaceEdgeId) -> FaceEdgeId>,
    )
}

impl<'a, V, E, F> Iterator for FaceRingEdges<'a, V, E, F> {
    type Item = (EdgeId, &'a Edge<E>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Single => None,

            Self::Multiple(iter) => {
                let plc = iter.common;
                iter.next().map(|(_, fe)| {
                    (fe.edge, &plc.edges[fe.edge])
                })
            }
        }
    }
}

impl<V, E, F> Index<VertexId> for Plc<V, E, F> {
    type Output = Vertex<V>;

    fn index(&self, index: VertexId) -> &Self::Output {
        self.vertex(index)
            .unwrap_or_else(|| panic!("Vertex {} does not exist", index))
    }
}

impl<V, E, F> Index<EdgeId> for Plc<V, E, F> {
    type Output = Edge<E>;

    fn index(&self, index: EdgeId) -> &Self::Output {
        self.edge(index)
            .unwrap_or_else(|| panic!("Edge {} does not exist", index))
    }
}

impl<V, E, F> Index<FaceId> for Plc<V, E, F> {
    type Output = Face<F>;

    fn index(&self, index: FaceId) -> &Self::Output {
        self.face(index)
            .unwrap_or_else(|| panic!("Face {} does not exist", index))
    }
}

/// Iterates over vertex ids and vertices.
pub type Vertices<'a, V> = id_map::Iter<'a, VertexId, Vertex<V>>;

/// Iterates over edge ids and edges.
pub type Edges<'a, E> = id_map::Iter<'a, EdgeId, Edge<E>>;

/// Iterates over face ids and faces.
pub type Faces<'a, F> = id_map::Iter<'a, FaceId, Face<F>>;

/// Iterates over vertex ids
pub type VertexIds<'a, V> = id_map::Keys<'a, VertexId, Vertex<V>>;

/// Keysates over edge ids
pub type EdgeIds<'a, E> = id_map::Keys<'a, EdgeId, Edge<E>>;

/// Keysates over face ids
pub type FaceIds<'a, F> = id_map::Keys<'a, FaceId, Face<F>>;

///// Iterates over vertex ids generated from extending the vertex list.
//pub type ExtendVertices<'a, V, I> = id_map::ExtendValues<'a, VertexId, Vertex<V>, I>;
//
///// Iterates over edge ids generated from extending the edge list.
//pub type ExtendEdges<'a, E, I> = id_map::ExtendValues<'a, EdgeId, Edge<E>, I>;
//
///// Iterates over face ids generated from extending the face list.
//pub type ExtendFaces<'a, F, I> = id_map::ExtendValues<'a, FaceId, Face<F>, I>;

#[cfg(test)]
mod tests {
    use super::*;
    use pathfinding::directed::bfs;

    /// Assert all the invariants of this PLC structure. This does not include
    /// not having intersections.
    #[cfg(test)]
    #[track_caller]
    fn assert_integrity<V, E, F>(plc: &Plc<V, E, F>) {
        // Vertex invariants
        for (id, vertex) in plc.vertices.iter() {
            if let Some(edge_out) = vertex.edge_out.valid() {
                assert_eq!(plc.edges[edge_out].vertices[0], id, "Out-edge {} of vertex {} does not point to the vertex.", edge_out, id);
            }

            if let Some(edge_in) = vertex.edge_in.valid() {
                assert_eq!(plc.edges[edge_in].vertices[1], id, "In-edge {} of vertex {} does not point to the vertex.", edge_in, id);
            }

            if let Some(ring) = vertex.ring.valid() {
                assert_eq!(plc.rings[ring].element, RingElement::Vertex(id),
                    "Face ring {} of vertex {} does not point to the vertex.", ring, id);
            }
        }

        // Edge invariants
        let mut edge_map = plc.edge_map.clone();
        for (id, edge) in plc.edges.iter() {
            assert!(edge.vertices[0].is_valid(), "Edge {} has an invalid source vertex", id);
            assert!(edge.vertices[1].is_valid(), "Edge {} has an invalid target vertex", id);
            assert_ne!(edge.vertices[0], edge.vertices[1], "Edge {} has vertices that are the same.", id);

            if let Some(edge_id) = edge_map.remove(&[edge.vertices[0], edge.vertices[1]]) {
                assert_eq!(edge_id, id, "Edge {} does not match id of edge {}-{} in edge map.", id, edge.vertices[0], edge.vertices[1]);
            } else {
                panic!("Edge {} ({}-{}) is missing in edge map.", id, edge.vertices[0], edge.vertices[1]);
            }

            if let Some(face_edge) = edge.face_edge.valid() {
                assert_eq!(plc.face_edges[face_edge].edge, id, "Face edge {} of edge {} does not point to the edge.", face_edge, id);
            }

            assert!(edge.next_edge_out.is_valid(), "Edge {} has an invalid next out-edge", id);
            assert!(edge.prev_edge_out.is_valid(), "Edge {} has an invalid prev out-edge", id);
            assert!(edge.next_edge_in.is_valid(), "Edge {} has an invalid next in-edge", id);
            assert!(edge.prev_edge_in.is_valid(), "Edge {} has an invalid prev in-edge", id);

            assert_eq!(plc.edges[edge.next_edge_out].vertices[0], edge.vertices[0],
                "Next out-edge {} of edge {} does not point to the same source vertex.", edge.next_edge_out, id);
            assert_eq!(plc.edges[edge.prev_edge_out].vertices[0], edge.vertices[0],
                "Prev out-edge {} of edge {} does not point to the same source vertex.", edge.prev_edge_out, id);
            assert_eq!(plc.edges[edge.next_edge_in].vertices[1], edge.vertices[1],
                "Next in-edge {} of edge {} does not point to the same target vertex.", edge.next_edge_in, id);
            assert_eq!(plc.edges[edge.prev_edge_in].vertices[1], edge.vertices[1],
                "Prev in-edge {} of edge {} does not point to the same target vertex.", edge.prev_edge_in, id);
        }
        assert_circular_lists(&plc.edges, |e| e.prev_edge_out, |e| e.next_edge_out, "out-edge");
        assert_circular_lists(&plc.edges, |e| e.prev_edge_in, |e| e.next_edge_in, "in-edge");

        if let Some(([v0, v1], edge)) = edge_map.iter().next() {
            panic!("Edge {} ({}-{}) is in edge map but is missing in edge list", edge, v0, v1);
        }

        // Face edge invariants
        for (id, face_edge) in plc.face_edges.iter() {
            assert!(face_edge.edge.is_valid(), "Face edge {} has an invalid edge", id);
            assert!(face_edge.ring.is_valid(), "Face edge {} has an invalid face ring", id);

            assert!(face_edge.next.is_valid(), "Face edge {} has an invalid next face edge", id);
            assert!(face_edge.prev.is_valid(), "Face edge {} has an invalid prev face edge", id);
            assert!(face_edge.next_on_edge.is_valid(), "Face edge {} has an invalid next face edge on the edge", id);
            assert!(face_edge.prev_on_edge.is_valid(), "Face edge {} has an invalid prev face edge on the edge", id);

            assert_eq!(plc.edges[plc.face_edges[face_edge.next].edge].vertices[0],
                plc.edges[face_edge.edge].vertices[1],
                "Next face edge {} of face edge {} does not start at the end of the face edge.", face_edge.next, id);
            assert_eq!(plc.edges[plc.face_edges[face_edge.prev].edge].vertices[1],
                plc.edges[face_edge.edge].vertices[0],
                "Prev face edge {} of face edge {} does not end at the start of the face edge.", face_edge.prev, id);

            assert_eq!(plc.face_edges[face_edge.next].ring, face_edge.ring,
                "Next face edge {} of face edge {} does not point to the same face ring.", face_edge.next, id);
            assert_eq!(plc.face_edges[face_edge.prev].ring, face_edge.ring,
                "Prev face edge {} of face edge {} does not point to the same face ring.", face_edge.prev, id);

            assert_eq!(plc.face_edges[face_edge.next_on_edge].edge, face_edge.edge,
                "Next face edge {} on the edge of face edge {} does not point to the same edge.", face_edge.next_on_edge, id);
            assert_eq!(plc.face_edges[face_edge.prev_on_edge].edge, face_edge.edge,
                "Prev face edge {} on the edge of face edge {} does not point to the same edge.", face_edge.prev_on_edge, id);
        }
        assert_circular_lists(&plc.face_edges, |fe| fe.prev, |fe| fe.next, "face edge");
        assert_circular_lists(&plc.face_edges, |fe| fe.prev_on_edge, |fe| fe.next_on_edge,
            "face edge (on edge)");

        // Face ring invariants
        for (id, ring) in plc.rings.iter() {
            match ring.element {
                RingElement::FaceEdge(face_edge) => {
                    assert!(face_edge.is_valid(), "Face ring {} has an invalid face edge", id);
                    assert_eq!(plc.face_edges[face_edge].ring, id, "Face edge {} of face ring {} does not point to the face ring.", face_edge, id);
                },

                RingElement::Vertex(vertex) => {
                    assert!(vertex.is_valid(), "Face ring {} has an invalid vertex", id);
                    assert_eq!(plc.vertices[vertex].ring, id, "Vertex {} of face ring {} does not point to the face ring.", vertex, id);
                },
            }
            assert!(ring.face.is_valid(), "Face ring {} has an invalid face.", id);

            assert!(ring.next.is_valid(), "Face ring {} has an invalid next face ring.", id);
            assert!(ring.prev.is_valid(), "Face ring {} has an invalid prev face ring.", id);

            assert_eq!(plc.rings[ring.next].face, ring.face, "Next face ring {} of face ring {} does not point to the same face.", ring.next, id);
            assert_eq!(plc.rings[ring.prev].face, ring.face, "Prev face ring {} of face ring {} does not point to the same face.", ring.prev, id);
        }
        assert_circular_lists(&plc.rings, |r| r.prev, |r| r.next, "face ring");

        // Face invariants
        for (id, face) in plc.faces.iter() {
            assert!(face.ring.is_valid(), "Face {} has an invalid face ring.", id);
            assert_eq!(plc.rings[face.ring].face, id, "Face ring {} of face {} does not point to the face.", face.ring, id);
        }
    }

    #[cfg(test)]
    #[track_caller]
    fn assert_circular_lists<Id: id_map::Id + Eq + Hash + Display + Debug, Val, PF: Fn(&Val) -> Id, NF: Fn(&Val) -> Id>
        (list: &IdMap<Id, Val>, prev_fn: PF, next_fn: NF, elem_name: &str)
    {
        let mut ids = list.keys().collect::<FnvHashSet<_>>();

        while let Some(id) = ids.iter().next().copied() {
            if let Some(mut bfs_loop) = bfs::bfs_loop(&id, |id| iter::once(next_fn(&list[*id]))) {
                for id in &bfs_loop {
                    ids.remove(id);
                }

                assert_eq!(bfs::bfs_loop(&id, |id| iter::once(prev_fn(&list[*id]))), {
                    bfs_loop.reverse();
                    Some(bfs_loop)
                }, "Reverse circular list of {0} {1} is not the reverse of the forward circular list.", elem_name, id);
            } else {
                panic!("Circular list of {0} {1} does not cycle back to the {0}.", elem_name, id);
            }
        }
    }

    #[track_caller]
    fn assert_vertices<V, E, F, I: IntoIterator<Item = VertexId>>(plc: &Plc<V, E, F>, vertices: I) {
        let result = plc.vertex_ids().collect::<FnvHashSet<_>>();
        let expect = vertices.into_iter().collect::<FnvHashSet<_>>();
        assert_eq!(result, expect, "Assert vertices failed");
    }

    #[track_caller]
    fn assert_edges<V, E, F, I: IntoIterator<Item = [VertexId; 2]>>(plc: &Plc<V, E, F>, edges: I) {
        let result = plc.edge_map.keys().copied().collect::<FnvHashSet<_>>();
        let expect = edges.into_iter().collect::<FnvHashSet<_>>();
        assert_eq!(result, expect, "Assert edges failed");
    }

    fn canonicalize_faces(faces: &mut [Vec<Vec<VertexId>>]) {
        for face in &mut *faces {
            for ring in &mut *face {
                // Canonicalize ring
                *ring = (0..ring.len())
                    .map(|i| {
                        let mut ring = ring.clone();
                        ring.rotate_left(i);
                        ring
                    }).min().unwrap();
            }

            // Canonicalize face
            face.sort();
        }

        // Canonicalize faces
        faces.sort();
    }

    #[track_caller]
    fn assert_faces<V, E, F, I: IntoIterator<Item = RI>, RI: IntoIterator<Item = VI>, VI: IntoIterator<Item = VertexId>>(plc: &Plc<V, E, F>, faces: I) {
        // This is tricky. The faces and rings could be in any order, but the vertices of each ring must be a cyclic permutation of each other.
        
        let mut result = plc.faces()
            .map(|(_, face)| face.rings(plc)
                .map(|(_, ring)| ring.vertices(plc).map(|(id, _)| id).collect::<Vec<_>>())
            .collect::<Vec<_>>()
        ).collect::<Vec<_>>();
        canonicalize_faces(&mut result);

        let mut expect = faces.into_iter()
            .map(|face| face.into_iter()
                .map(|ring| ring.into_iter().collect::<Vec<_>>()
            ).collect::<Vec<_>>()
        ).collect::<Vec<_>>();
        canonicalize_faces(&mut expect);

        assert_eq!(result, expect, "Assert faces failed");
    }

    fn v_ids(ids: Vec<u32>) -> Vec<VertexId> {
        ids.into_iter().map(|id| VertexId(id)).collect()
    }

    fn e_ids(ids: Vec<[u32; 2]>) -> Vec<[VertexId; 2]> {
        ids.into_iter().map(|[a, b]| [VertexId(a), VertexId(b)]).collect()
    }

    fn f_id(ids: Vec<Vec<u32>>) -> Vec<Vec<VertexId>> {
        ids.into_iter().map(|ids| v_ids(ids)).collect()
    }

    fn f_ids(ids: Vec<Vec<Vec<u32>>>) -> Vec<Vec<Vec<VertexId>>> {
        ids.into_iter().map(|ids| ids.into_iter().map(|ids| v_ids(ids)).collect()).collect()
    }

    fn plc_from_vef(num_vertices: usize, edges: Vec<[u32; 2]>, faces: Vec<Vec<Vec<u32>>>) -> Plc<(), (), ()> {
        let mut plc = Plc::new(|| (), || (), || ());
        (0..num_vertices)
            .map(|_| plc.add_vertex(Pt3::origin(), ())).for_each(|_| {});

        e_ids(edges).into_iter().map(|vertices| plc.add_edge(vertices, ())).for_each(|_| {});
        f_ids(faces).into_iter().map(|face| plc.add_face(face, ())).for_each(|_| {});

        plc
    }

    fn plc_from_pos_vef(points: Vec<Pt3>, edges: Vec<[u32; 2]>, faces: Vec<Vec<Vec<u32>>>) -> Plc<(), (), ()> {
        let mut plc = Plc::new(|| (), || (), || ());
        points.into_iter().map(|point| plc.add_vertex(point, ())).for_each(|_| {});

        e_ids(edges).into_iter().map(|vertices| plc.add_edge(vertices, ())).for_each(|_| {});
        f_ids(faces).into_iter().map(|face| plc.add_face(face, ())).for_each(|_| {});

        plc
    }

    #[test]
    fn test_empty() {
        let plc = Plc::new(|| (), || (), || ());
        assert_integrity(&plc);
        assert_vertices(&plc, vec![]);
        assert_edges(&plc, vec![]);
        assert_faces(&plc, vec![] as Vec<Vec<Vec<VertexId>>>);
    }

    #[test]
    fn test_add_vertex() {
        let mut plc = Plc::new(|| 0, || (), || ());
        let id = plc.add_vertex(Pt3::new(1.0, 2.0, 3.0), 1);

        assert_integrity(&plc);
        assert_vertices(&plc, vec![id]);
        assert_edges(&plc, vec![]);
        assert_faces(&plc, vec![] as Vec<Vec<Vec<VertexId>>>);
        assert_eq!(plc[id].position(), Pt3::new(1.0, 2.0, 3.0));
        assert_eq!(*plc[id].value(), 1);
    }

    #[test]
    fn test_add_edge() {
        let mut plc = Plc::new(|| (), || 0, || ());
        let v = (0..3)
            .map(|_| plc.add_vertex(Pt3::origin(), ()))
            .collect::<Vec<_>>();

        let id = plc.add_edge([v[0], v[1]], 1);
        
        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2]));
        assert_edges(&plc, e_ids(vec![[0, 1]]));
        assert_faces(&plc, vec![] as Vec<Vec<Vec<VertexId>>>);
        assert_eq!(*plc[id].value(), 1);

        let id = plc.add_edge_if_absent([v[0], v[1]], 2);
        assert_eq!(*plc[id].value(), 1);

        plc.add_edge([v[1], v[2]], 1);
        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2]]));
        assert_faces(&plc, vec![] as Vec<Vec<Vec<VertexId>>>);
    }

    #[test]
    fn test_add_edge_all() {
        let mut plc = Plc::new(|| (), || (), || ());
        let v = (0..3)
            .map(|_| plc.add_vertex(Pt3::origin(), ()))
            .collect::<Vec<_>>();

        e_ids(vec![
            [0, 1], [1, 0], [1, 2], [2, 1], [2, 0], [0, 2]
        ]).into_iter().map(|vertices| plc.add_edge(vertices, ())).for_each(|_| {});

        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2], [2, 0], [1, 0], [2, 1], [0, 2]]));
        assert_faces(&plc, vec![] as Vec<Vec<Vec<VertexId>>>);
    }

    #[test]
    fn test_add_triangle() {
        let mut plc = Plc::new(|| (), || (), || 0);
        let v = (0..4)
            .map(|_| plc.add_vertex(Pt3::origin(), ()))
            .collect::<Vec<_>>();

        // Triangle
        let id = plc.add_face(f_id(vec![vec![1, 2, 3]]), 3);
        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2, 3]));
        assert_edges(&plc, e_ids(vec![[3, 1], [1, 2], [2, 3]]));
        assert_faces(&plc, f_ids(vec![vec![vec![2, 3, 1]]]));
        assert_eq!(*plc[id].value(), 3);
    }

    #[test]
    fn test_add_mwb_polyhedron() {
        let plc = plc_from_vef(6, vec![], vec![
            vec![vec![0, 1, 2]],
            vec![vec![0, 2, 4, 3]],
            vec![vec![5, 4, 2, 1]],
        ]);

        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2, 3, 4, 5]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2], [2, 0], [0, 2], [2, 4], [4, 3], [3, 0], [5, 4], [4, 2], [2, 1], [1, 5]]));
        assert_faces(&plc, f_ids(vec![
            vec![vec![1, 2, 0]],
            vec![vec![2, 1, 5, 4]],
            vec![vec![0, 2, 4, 3]],
        ]));
    }

    #[test]
    fn test_add_non_manifold_faces() {
        let plc = plc_from_vef(4, vec![], vec![
            vec![vec![0, 1, 2]],
            vec![vec![0, 1, 3]],
        ]);

        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2, 3]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2], [2, 0], [1, 3], [3, 0]]));
        assert_faces(&plc, f_ids(vec![
            vec![vec![0, 1, 2]],
            vec![vec![0, 1, 3]],
        ]));
    }

    #[test]
    fn test_add_face_with_holes() {
        let plc = plc_from_vef(8, vec![], vec![
            vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]],
            vec![vec![7, 6, 5, 4]],
        ]);

        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2, 3, 4, 5, 6, 7]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [7, 6], [6, 5], [5, 4], [4, 7]]));
        assert_faces(&plc, f_ids(vec![
            vec![vec![6, 5, 4, 7]],
            vec![vec![6, 7, 4, 5], vec![0, 1, 2, 3]],
        ]));
    }

    #[test]
    fn test_add_face_with_isolated_vertex() {
        let plc = plc_from_vef(8, vec![], vec![
            vec![vec![0, 1, 2, 3], vec![4]],
        ]);

        assert_integrity(&plc);
        assert_vertices(&plc, v_ids(vec![0, 1, 2, 3, 4, 5, 6, 7]));
        assert_edges(&plc, e_ids(vec![[0, 1], [1, 2], [2, 3], [3, 0]]));
        assert_faces(&plc, f_ids(vec![
            vec![vec![0, 1, 2, 3], vec![4]],
        ]));
    }

    #[test]
    #[cfg(feature = "obj")]
    fn test_load_obj() {
        let plc = Plc::load_obj("assets/monkey.obj", || (), || (), || ()).unwrap();
        assert_integrity(&plc);
    }

    #[test]
    fn test_invalid_degenerate_face() {
        let plc = plc_from_pos_vef(vec![
            Pt3::new(0.0, 0.0, 0.0)
        ], vec![], vec![vec![vec![0]]]);

        assert_eq!(plc.validate(0.1).unwrap_err().reason(), &ValidateReason::DegenerateFace(vec![FaceId(0)]));
    }

    #[test]
    fn test_invalid_twisted() {
        let plc = plc_from_pos_vef(vec![
            Pt3::new(0.0, 0.0, 0.1),
            Pt3::new(0.0, 1.0, 0.0),
            Pt3::new(1.0, 1.0, 0.1),
            Pt3::new(1.0, 0.0, 0.0),
        ], vec![], vec![vec![vec![0, 1, 2, 3]]]);

        // Note: Coplanarity tolerance is half distance tolerance.
        plc.validate(0.06).unwrap_err();
    }

    #[test]
    fn test_invalid_intersects() {
        let plc = plc_from_pos_vef(vec![
            Pt3::new(0.0, 0.0, 0.0),
            Pt3::new(0.0, 0.0, 3.0),
            Pt3::new(0.0, 3.0, 3.0),
            Pt3::new(0.0, 3.0, 0.0),
            Pt3::new(0.0, 1.0, 1.0),
            Pt3::new(0.0, 1.0, 2.0),
            Pt3::new(0.0, 2.0, 2.0),
            Pt3::new(0.0, 2.0, 1.0),
            Pt3::new(0.0, 0.5, 1.0),
        ], vec![], vec![vec![vec![0, 1, 2, 3], vec![7, 6, 5, 4]]]);

        assert_eq!(plc.validate(0.1).unwrap_err().reason(),
            &ValidateReason::Intersects(vec![(Element::Vertex(VertexId(8)), Element::Face(FaceId(0)), 0.0)]));
    }

    #[test]
    fn test_valid_isolated_vertex() {
        let plc = plc_from_pos_vef(vec![
            Pt3::new(0.0, 0.0, 0.0),
            Pt3::new(0.0, 0.0, 3.0),
            Pt3::new(0.0, 3.0, 3.0),
            Pt3::new(0.0, 3.0, 0.0),
            Pt3::new(0.0, 1.0, 1.0),
            Pt3::new(0.0, 1.0, 2.0),
            Pt3::new(0.0, 2.0, 2.0),
            Pt3::new(0.0, 2.0, 1.0),
            Pt3::new(0.0, 0.5, 1.0),
        ], vec![], vec![vec![vec![0, 1, 2, 3], vec![7, 6, 5, 4], vec![8]]]);

        assert!(plc.validate(0.1).is_ok());
    }
}