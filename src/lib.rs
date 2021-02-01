mod id_map;
pub mod tet;

use nalgebra::Point3;

type Pt3 = Point3<f64>;

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
        }
    };
}

id!(pub struct VertexId);
id!(pub struct TetId);

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
        /// This is an alias.
        $pub fn $alias$(<$($gen),*>)?($($arg: $arg_type),*) -> $ret {
            Self::$name $(::<$($gen),*>)?($($arg),*)
        }
    }
}
