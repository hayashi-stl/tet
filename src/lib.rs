mod id_map;
pub mod tet;

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
            pub(crate) fn invalid() -> Self {
                Self(crate::id_map::IdType::MAX)
            }
        }
    };
}

id!(pub struct VertexId);
id!(pub struct TetId);