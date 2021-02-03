use fnv::FnvHashSet;
use std::iter::Enumerate;
use std::vec;
use std::{
    iter::FromIterator,
    marker::PhantomData,
    ops::{Index, IndexMut},
    slice,
};

pub type IdType = u32;

pub trait Id: Copy {
    fn int(self) -> IdType;

    fn from_int(int: IdType) -> Self;
}

/// Map that reuses keys and allows specific key insertion.
/// Expects few deletions because they would make iteration slower.
#[derive(Clone, Debug)]
pub(crate) struct IdMap<K, V> {
    map: Vec<Option<V>>,
    free: FnvHashSet<IdType>,
    marker: PhantomData<K>,
}

impl<K, V> Default for IdMap<K, V> {
    fn default() -> Self {
        Self {
            map: vec![],
            free: FnvHashSet::default(),
            marker: PhantomData,
        }
    }
}

impl<K: Id, V> IdMap<K, V> {
    /// Constructs a new id map.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Gets the length of the map
    pub(crate) fn len(&self) -> usize {
        self.map.len() - self.free.len()
    }

    /// Gets the value with a specific key, if it exists.
    pub(crate) fn get(&self, key: K) -> Option<&V> {
        self.map.get(key.int() as usize).and_then(|v| v.as_ref())
    }

    /// Gets the value with a specific key mutably, if it exists.
    pub(crate) fn get_mut(&mut self, key: K) -> Option<&mut V> {
        self.map
            .get_mut(key.int() as usize)
            .and_then(|v| v.as_mut())
    }

    /// Inserts a value into the map and returns the key for that value.
    pub(crate) fn insert(&mut self, value: V) -> K {
        if self.free.is_empty() {
            self.map.push(Some(value));
            K::from_int(self.map.len() as IdType - 1)
        } else {
            let key = self.free.iter().next().copied().unwrap();
            self.free.remove(&key);
            self.map[key as usize] = Some(value);
            K::from_int(key)
        }
    }

    /// Inserts a value into the map with a specific key and returns the value that was previously there, if it existed.
    pub(crate) fn insert_with_key(&mut self, key: K, value: V) -> Option<V> {
        if key.int() as usize >= self.map.len() {
            // Create more space
            self.free.extend(self.map.len() as IdType..key.int());
            self.map.resize_with(key.int() as usize + 1, || None);
            self.map[key.int() as usize] = Some(value);
            None
        } else {
            let old = std::mem::replace(&mut self.map[key.int() as usize], Some(value));
            if old.is_none() {
                self.free.remove(&key.int());
            }
            old
        }
    }

    /// Extends the values with an iterator and returns an iterator over the keys in order of insertion.
    #[must_use = "This is an iterator."]
    pub(crate) fn extend_values<I: IntoIterator<Item = V>>(
        &mut self,
        values: I,
    ) -> ExtendValues<K, V, I::IntoIter> {
        ExtendValues {
            map: self,
            values: values.into_iter(),
        }
    }

    /// Removes the value with a specific key and returns that value, if it existed.
    pub(crate) fn remove(&mut self, key: K) -> Option<V> {
        let opt_value = self
            .map
            .get_mut(key.int() as usize)
            .unwrap_or(&mut None)
            .take();

        if opt_value.is_some() {
            self.free.insert(key.int());
        }

        opt_value
    }

    /// Keeps only the values that satisfy some predicate.
    pub(crate) fn retain<F: FnMut(K, &V) -> bool>(&mut self, mut pred: F) {
        for i in 0..self.map.len() {
            let key = K::from_int(i as IdType);
            if !self.get(key).map(|v| pred(key, v)).unwrap_or(true) {
                self.remove(key);
            }
        }
    }

    /// Clears the map.
    pub(crate) fn clear(&mut self) {
        self.map.clear();
        self.free.clear();
    }

    /// Iterates over (key, &value) pairs.
    pub(crate) fn iter(&self) -> Iter<K, V> {
        Iter {
            iter: self.map.iter().enumerate(),
            marker: PhantomData,
        }
    }

    /// Iterates over (key, &mut value) pairs.
    pub(crate) fn iter_mut(&mut self) -> IterMut<K, V> {
        IterMut {
            iter: self.map.iter_mut().enumerate(),
            marker: PhantomData,
        }
    }

    /// Iterates over keys.
    pub(crate) fn keys(&self) -> Keys<K, V> {
        Keys {
            iter: self.map.iter().enumerate(),
            marker: PhantomData,
        }
    }

    /// Turns this into an iterator over the values.
    pub(crate) fn into_values(self) -> IntoValues<V> {
        IntoValues {
            iter: self.map.into_iter(),
        }
    }

    /// Iterates over values.
    pub(crate) fn values(&self) -> Values<V> {
        Values {
            iter: self.map.iter(),
        }
    }

    /// Iterates over values mutably.
    pub(crate) fn values_mut(&mut self) -> ValuesMut<V> {
        ValuesMut {
            iter: self.map.iter_mut(),
        }
    }
}

impl<K: Id, V> IntoIterator for IdMap<K, V> {
    type Item = (K, V);
    type IntoIter = IntoIter<K, V>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            iter: self.map.into_iter().enumerate(),
            marker: PhantomData,
        }
    }
}

impl<K: Id, V> Extend<(K, V)> for IdMap<K, V> {
    /// Extend with (key, value) pairs
    fn extend<T: IntoIterator<Item = (K, V)>>(&mut self, iter: T) {
        for (k, v) in iter.into_iter() {
            self.insert_with_key(k, v);
        }
    }
}

impl<K: Id, V> FromIterator<(K, V)> for IdMap<K, V> {
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        let mut map = Self::new();
        map.extend(iter);
        map
    }
}

impl<K: Id, V> FromIterator<V> for IdMap<K, V> {
    /// The keys are the ones corresponding to indexes 0..n, where n is the length of the iterator.
    /// The keys are in order of the values returned by the iterator.
    fn from_iter<T: IntoIterator<Item = V>>(iter: T) -> Self {
        let mut map = Self::new();
        map.map.extend(iter.into_iter().map(|v| Some(v)));
        map
    }
}

impl<K: Id, V> Index<K> for IdMap<K, V> {
    type Output = V;

    fn index(&self, index: K) -> &Self::Output {
        self.get(index)
            .unwrap_or_else(|| panic!("Index {} does not exist.", index.int()))
    }
}

impl<K: Id, V> IndexMut<K> for IdMap<K, V> {
    fn index_mut(&mut self, index: K) -> &mut Self::Output {
        self.get_mut(index)
            .unwrap_or_else(|| panic!("Index {} does not exist.", index.int()))
    }
}

pub struct ExtendValues<'a, K, V, I> {
    map: &'a mut IdMap<K, V>,
    values: I,
}

impl<'a, K: Id, V, I: Iterator<Item = V>> Iterator for ExtendValues<'a, K, V, I> {
    type Item = K;

    fn next(&mut self) -> Option<Self::Item> {
        self.values.next().map(|value| self.map.insert(value))
    }
}

pub(crate) struct Keys<'a, K, V> {
    iter: Enumerate<slice::Iter<'a, Option<V>>>,
    marker: PhantomData<K>,
}

impl<'a, K: Id, V> Iterator for Keys<'a, K, V> {
    type Item = K;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((key, value)) = self.iter.next() {
            if let Some(_) = value.as_ref() {
                return Some(K::from_int(key as IdType));
            }
        }
        None
    }
}

pub(crate) struct IntoValues<V> {
    iter: vec::IntoIter<Option<V>>,
}

impl<V> Iterator for IntoValues<V> {
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(value) = self.iter.next() {
            if let Some(value) = value {
                return Some(value);
            }
        }
        None
    }
}

pub(crate) struct Values<'a, V> {
    iter: slice::Iter<'a, Option<V>>,
}

impl<'a, V> Iterator for Values<'a, V> {
    type Item = &'a V;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(value) = self.iter.next() {
            if let Some(value) = value.as_ref() {
                return Some(value);
            }
        }
        None
    }
}

pub(crate) struct ValuesMut<'a, V> {
    iter: slice::IterMut<'a, Option<V>>,
}

impl<'a, V> Iterator for ValuesMut<'a, V> {
    type Item = &'a mut V;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(value) = self.iter.next() {
            if let Some(value) = value.as_mut() {
                return Some(value);
            }
        }
        None
    }
}

pub(crate) struct IntoIter<K, V> {
    iter: Enumerate<vec::IntoIter<Option<V>>>,
    marker: PhantomData<K>,
}

impl<K: Id, V> Iterator for IntoIter<K, V> {
    type Item = (K, V);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, value)) = self.iter.next() {
            if let Some(value) = value {
                return Some((K::from_int(i as IdType), value));
            }
        }
        None
    }
}

pub struct Iter<'a, K, V> {
    iter: Enumerate<slice::Iter<'a, Option<V>>>,
    marker: PhantomData<K>,
}

impl<'a, K: Id, V> Iterator for Iter<'a, K, V> {
    type Item = (K, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, value)) = self.iter.next() {
            if let Some(value) = value.as_ref() {
                return Some((K::from_int(i as IdType), value));
            }
        }
        None
    }
}

pub(crate) struct IterMut<'a, K, V> {
    iter: Enumerate<slice::IterMut<'a, Option<V>>>,
    marker: PhantomData<K>,
}

impl<'a, K: Id, V> Iterator for IterMut<'a, K, V> {
    type Item = (K, &'a mut V);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, value)) = self.iter.next() {
            if let Some(value) = value.as_mut() {
                return Some((K::from_int(i as IdType), value));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl Id for u32 {
        fn int(self) -> IdType {
            self
        }

        fn from_int(int: IdType) -> Self {
            int
        }
    }

    #[test]
    fn test_default() {
        let map = IdMap::<u32, u32>::default();
        assert_eq!(map.len(), 0);
    }

    #[test]
    fn test_insert_one() {
        let mut map = IdMap::<u32, _>::new();
        let key = map.insert(1);
        assert_eq!(map[key], 1);
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn test_insert_multiple() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        let key2 = map.insert(6);
        assert_eq!(map[key0], 2);
        assert_eq!(map[key1], 5);
        assert_eq!(map[key2], 6);
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn test_remove_out_of_bounds() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        assert_eq!(map.remove(2), None);
        assert_eq!(map[key0], 2);
        assert_eq!(map[key1], 5);
        assert_eq!(map.len(), 2);
    }

    #[test]
    fn test_remove_present() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        assert_eq!(map.remove(key0), Some(2));
        assert_eq!(map[key1], 5);
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn test_remove_absent() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        assert_eq!(map.remove(key0), Some(2));
        assert_eq!(map.remove(key0), None);
        assert_eq!(map[key1], 5);
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn test_insert_with_key_resize() {
        let mut map = IdMap::<u32, _>::new();
        let key = map.insert(1);
        assert_eq!(map.insert_with_key(1, 8), None);
        assert_eq!(map.insert_with_key(4, 9), None);
        assert_eq!(map[key], 1);
        assert_eq!(map[1], 8);
        assert_eq!(map[4], 9);
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn test_insert_with_key_replace() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        assert_eq!(map.insert_with_key(key0, 6), Some(2));
        assert_eq!(map[key0], 6);
        assert_eq!(map[key1], 5);
        assert_eq!(map.len(), 2);
    }

    #[test]
    fn test_insert_with_key_fill_hole() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        assert_eq!(map.remove(key0), Some(2));
        assert_eq!(map.insert_with_key(key0, 6), None);
        assert_eq!(map[key0], 6);
        assert_eq!(map[key1], 5);
        assert_eq!(map.len(), 2);
    }

    #[test]
    fn test_into_values() {
        let mut map = IdMap::<u32, _>::new();
        let _ey0 = map.insert(2);
        let key1 = map.insert(5);
        let _ey2 = map.insert(6);
        map.remove(key1);
        let mut iter = map.into_values();
        assert_eq!(iter.next(), Some(2));
        assert_eq!(iter.next(), Some(6));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_into_iterator() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        let key2 = map.insert(6);
        map.remove(key1);
        let mut iter = map.into_iter();
        assert_eq!(iter.next(), Some((key0, 2)));
        assert_eq!(iter.next(), Some((key2, 6)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_extend() {
        let mut map = IdMap::<u32, _>::new();
        let key0 = map.insert(2);
        let key1 = map.insert(5);
        let _ey2 = map.insert(6);
        map.remove(key1);
        map.extend(vec![(1, 1), (2, 4), (3, 9), (5, 3)]);
        assert_eq!(map[key0], 2);
        assert_eq!(map[1], 1);
        assert_eq!(map[2], 4);
        assert_eq!(map[3], 9);
        assert_eq!(map[5], 3);
        assert_eq!(map.len(), 5);
    }

    #[test]
    fn test_extend_values() {
        let mut map = IdMap::<u32, _>::new();
        let _ey0 = map.insert(2);
        let key1 = map.insert(5);
        let _ey2 = map.insert(6);
        map.remove(key1);
        let keys = map.extend_values(vec![1, 4, 9, 3]);
        assert_eq!(keys.collect::<Vec<_>>(), vec![1, 3, 4, 5]);
    }
}
