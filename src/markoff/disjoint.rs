use either::{Either, Left, Right};
use std::collections::{HashMap, HashSet};

/// A set of disjoint sets of instances of `K`.
pub struct Disjoint<K> {
    disjoint: HashMap<K, Either<K, u128>>,
    orbits: HashSet<K>,
}

impl<K: Eq + Clone + std::hash::Hash> Disjoint<K> {
    /// Creates a new, empty set of disjoint sets.
    pub fn new() -> Disjoint<K> {
        Disjoint {
            disjoint: HashMap::new(),
            orbits: HashSet::new(),
        }
    }

    /// Returns an `Iterator` yielding, for each disjoint set, a representative and the size of the
    /// disjoint set, respectively.
    /// There are no guarantees about which element will be the chosen representative.
    pub fn get_sets(&self) -> impl Iterator<Item = (&K, u128)> {
        self.orbits
            .iter()
            .map(|key| self.disjoint.get(key).map(|e| (key, e)).unwrap())
            .map(|(k, e)| e.as_ref().right().map(|d| (k, *d)).unwrap())
    }

    /// Merge two disjoint sets; specifically, if `one` is in $S$ and `two` is in $T$, then
    /// $S$ and $T$ are replaced by $S \cup T$.
    /// If either `one` or `two` were not previously in this set of disjoint sets, then they are
    /// considered to be singletons, and then merged.
    pub fn associate(&mut self, one: K, two: K) {
        match (self.root(&one), self.root(&two)) {
            (None, None) => {
                let count = if one == two { 1 } else { 2 };

                self.orbits.insert(one.clone());
                self.disjoint.insert(two, Left(one.clone()));
                self.disjoint.insert(one, Right(count));
            }
            (Some((k, d)), None) => {
                self.disjoint.insert(two, Left(k.clone()));
                self.disjoint.insert(k, Right(d + 1));
            }
            (None, Some((k, d))) => {
                self.disjoint.insert(one, Left(k.clone()));
                self.disjoint.insert(k, Right(d + 1));
            }
            (Some((k1, d1)), Some((k2, d2))) if k1 != k2 => {
                if d1 >= d2 {
                    self.disjoint.insert(k2.clone(), Left(k1.clone()));
                    self.disjoint.insert(k1, Right(d1 + d2));
                    self.orbits.remove(&k2);
                } else {
                    self.disjoint.insert(k1.clone(), Left(k2.clone()));
                    self.disjoint.insert(k2, Right(d1 + d2));
                    self.orbits.remove(&k1);
                }
            },
            _ => {}
        }
    }

    fn root(&self, key: &K) -> Option<(K, u128)> {
        match self.disjoint.get(key) {
            None => None,
            Some(Left(k2)) => self.root(k2),
            Some(Right(d)) => Some((key.clone(), *d)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assoc() {
        let mut disjoint: Disjoint<u32> = Disjoint::new();
        let assocs = vec![(1, 2), (2, 3), (4, 5), (6, 7), (8, 9), (6, 2), (9, 4)];
        for (x, y) in assocs {
            disjoint.associate(x, y);
        }
        let orbits: Vec<(&u32, u128)> = disjoint.get_sets().collect();
        assert_eq!(orbits.len(), 2);
    }
}
