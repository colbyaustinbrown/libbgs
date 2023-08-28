use either::{Either, Left, Right};
use std::collections::{HashMap, HashSet};

pub struct Disjoint<K> {
    disjoint: HashMap<K, Either<K, u128>>,
    orbits: HashSet<K>,
}

impl<K: Eq + Clone + std::hash::Hash> Disjoint<K> {
    pub fn new() -> Disjoint<K> {
        Disjoint {
            disjoint: HashMap::new(),
            orbits: HashSet::new(),
        }
    }

    pub fn get_sets(&self) -> impl Iterator<Item = (&K, u128)> {
        self.orbits
            .iter()
            .filter_map(|key| self.disjoint.get(key).map(|e| (key, e)))
            .filter_map(|(k, e)| e.as_ref().right().map(|d| (k, *d)))
    }

    pub fn associate(&mut self, one: K, two: K) {
        match (self.root(&one), self.root(&two)) {
            (None, None) => {
                self.orbits.insert(one.clone());
                self.disjoint.insert(two, Left(one.clone()));
                self.disjoint.insert(one, Right(2));
            }
            (Some((k, d)), None) => {
                self.disjoint.insert(two, Left(k.clone()));
                self.disjoint.insert(one.clone(), Right(d + 1));
            }
            (None, Some((k, d))) => {
                self.disjoint.insert(one, Left(k.clone()));
                self.disjoint.insert(two.clone(), Right(d + 1));
            }
            (Some((k1, d1)), Some((k2, d2))) => {
                if d1 >= d2 {
                    let k1_ = k1.clone();
                    self.disjoint.insert(k2.clone(), Left(k1.clone()));
                    self.disjoint.insert(k1_, Right(d1 + d2));
                } else {
                    let k2_ = k2.clone();
                    self.disjoint.insert(k1.clone(), Left(k2.clone()));
                    self.disjoint.insert(k2_, Right(d1 + d2));
                }
            }
        }
    }

    fn root<'a>(&'a self, key: &'a K) -> Option<(&'a K, u128)> {
        match self.disjoint.get(key) {
            None => None,
            Some(Left(k2)) => self.root(k2),
            Some(Right(d)) => Some((key, *d)),
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
