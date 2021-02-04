mod id_map;
pub mod plc;
pub mod tet;
mod util;

pub use crate::plc::Plc;
pub use crate::tet::{TetMesh, TetWalker};

use nalgebra::{Point1, Point3, Vector1, Vector3};

type Pt1 = Point1<f64>;
type Pt3 = Point3<f64>;
type Vec1 = Vector1<f64>;
type Vec3 = Vector3<f64>;

#[macro_export]
#[doc(hidden)]
macro_rules! id {
    ($(#[$attr:meta])* $pub:vis struct $name:ident) => {
        $(#[$attr])*
        #[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
        $pub struct $name(pub(crate) crate::id_map::IdType);

        impl crate::id_map::Id for $name {
            fn int(self) -> crate::id_map::IdType {
                self.0
            }

            fn from_int(int: crate::id_map::IdType) -> Self {
                Self(int)
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                std::fmt::Display::fmt(&self.0, f)
            }
        }

        impl $name {
            /// Returns an invalid id.
            pub(crate) const fn invalid() -> Self {
                Self(crate::id_map::IdType::MAX)
            }

            /// Returns this id if it's valid or None otherwise.
            #[allow(dead_code)]
            pub(crate) fn valid(self) -> Option<Self> {
                if self == Self::invalid() {
                    None
                } else {
                    Some(self)
                }
            }

            /// Returns whether this id is valid.
            #[allow(dead_code)]
            pub(crate) fn is_valid(self) -> bool {
                self != Self::invalid()
            }
        }
    };
}

id! {
    /// A vertex id, both for tet meshes and PLCs
    pub struct VertexId
}

#[macro_export]
#[doc(hidden)]
macro_rules! alias {
    (
        $alias:ident,
        $(#[$attr:meta])*
        $pub:vis fn $name:ident$(<$($gen:ident),*>)?($($(($mut:tt))? $arg:ident: $arg_type:ty),*) -> $ret:ty {
            $($body:tt)*
        }
    ) => {
        $(#[$attr])*
        $pub fn $name$(<$($gen),*>)?($($($mut)? $arg: $arg_type),*) -> $ret {
            $($body)*
        }

        $(#[$attr])*
        #[inline(always)]
        /// This is an alias.
        $pub fn $alias$(<$($gen),*>)?($($arg: $arg_type),*) -> $ret {
            Self::$name $(::<$($gen),*>)?($($arg),*)
        }
    }
}
