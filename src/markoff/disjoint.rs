use std::collections::{HashSet, HashMap};
use std::rc::Rc;
use std::cmp::Ordering;
use either::{Left, Right, Either};
use std::cell::{RefCell};

#[derive(Debug)]
struct Set<V> {
    rank: u16,
    data: V
}

#[derive(Debug)]
struct SetPtr<V> {
    ptr: Rc<RefCell<Either<SetPtr<V>, Set<V>>>>
}

#[derive(Debug)]
struct Disjoint<K, V> {
    disjoint: HashMap<K, SetPtr<V>>,
    orbits: HashSet<K>
}

impl<K: Eq + Clone + std::hash::Hash, V> Disjoint<K, V> {
    pub fn new() -> Disjoint<K, V> {
        Disjoint {
            disjoint: HashMap::new(),
            orbits: HashSet::new()
        }
    }

    pub fn get_orbits<'a>(&'a self) -> impl Iterator<Item = &'a Set<V>> {
        self.orbits.iter()
            .filter_map(|key| {
                self.disjoint.get(key)
            })
            .filter_map(|e| {
                let ptr = unsafe { (&*e.ptr.as_ptr()).as_ref() };
                ptr.right()
            })
    }

    pub fn associate_with(&mut self, one: &K, two: &K, data: V) {
        match (self.disjoint.get(one), self.disjoint.get(two)) {
            (None, None) => {
                let op = SetPtr::new(Set::new(data));
                self.orbits.insert(one.clone());
                self.disjoint.insert(two.clone(), op.point());
                self.disjoint.insert(one.clone(), op);
            },
            (Some(x), None) => {
                self.disjoint.insert(two.clone(), x.point());
            },
            (None, Some(y)) => {
                self.disjoint.insert(one.clone(), y.point());
            },
            (Some(x), Some(y)) => {
                let (p1,o1) = x.root();
                let (p2,o2) = y.root();
                match o1.rank.cmp(&o2.rank) {
                    Ordering::Less => {
                        p1.ptr.replace(Left(p2.point()));
                    },
                    Ordering::Greater => {
                        p2.ptr.replace(Left(p1.point()));
                    },
                    Ordering::Equal => {
                        self.orbits.remove(one);
                        p1.ptr.replace(Left(p2.point()));
                        p2.ptr.replace(Right(Set {
                            rank: o2.rank,
                            data
                        }));
                    }
                }
            }
        }
    }
}

impl<K: Eq + Clone + std::hash::Hash, V: Default> Disjoint<K, V> {
    fn associate(&mut self, one: &K, two: &K) {
        self.associate_with(one, two, V::default());
    }
}

impl<V> Set<V> {
    fn new(data: V) -> Set<V> {
        Set {
            rank: 0,
            data
        }
    }
}

impl<V: Default> Default for Set<V> {
    fn default() -> Set<V> {
        Set {
            rank: 0,
            data: V::default()
        }
    }
}

impl<V> SetPtr<V> {
    fn new(orbit: Set<V>) -> SetPtr<V> {
        SetPtr { 
            ptr: Rc::new(RefCell::new(Right(orbit)))
        }
    }

    fn point(&self) -> SetPtr<V> {
        SetPtr {
            ptr: Rc::clone(&self.ptr)
        }
    }

    fn root(&self) -> (&SetPtr<V>, &Set<V>) {
        unsafe {
            match &*self.ptr.as_ptr() {
                Left(l) => {
                    l.root()
                },
                Right(r) => {
                    (self, r)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assoc() {
        let mut disjoint: Disjoint<u32, ()> = Disjoint::new(); 
        let assocs = vec![(1, 2), (2, 3), (4, 5), (6, 7), (8, 9), (6, 2), (9, 4)];
        for (x,y) in assocs {
            disjoint.associate_with(&x, &y, ());
        }
        let orbits: Vec<&Set<()>> = disjoint.get_orbits().collect();
        assert_eq!(orbits.len(), 2);
    }
}
