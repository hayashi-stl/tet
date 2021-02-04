//! Loads a PLC, triangulates its vertices, and exports the result as an obj
//! where each tet has its own vertices.
//! Requires the "obj" feature.

use std::env;
use tet::{Plc, TetMesh};

fn print_usage() {
    println!("Usage: <program> <input_obj> <output_obj>");
}

fn main() {
    #[cfg(feature = "obj")]
    {
        let mut args = env::args();
        args.next();
        let input = args.next().unwrap_or_else(|| {
            print_usage();
            panic!("Missing input obj")
        });
        let output = args.next().unwrap_or_else(|| {
            print_usage();
            panic!("Missing output obj")
        });

        let plc = Plc::load_obj(input, || (), || (), || ()).expect("Could not load input");
        if let Err(err) = plc.validate(0.0015) {
            panic!("{}", err)
        }

        let tets = TetMesh::delaunay_from_vertices(
            plc.vertices().map(|(_, v)| (v.position(), *v.value())),
            || ()
        );

        tets.export_debug_obj(output).expect("Coult not save output");
    }
    #[cfg(not(feature = "obj"))]
    panic!("This example requires the \"obj\" feature.");
}