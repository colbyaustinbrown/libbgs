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
struct Forest<K, V> {
    forest: HashMap<K, SetPtr<V>>,
    orbits: HashSet<K>
}

impl<K: Eq + Clone + std::hash::Hash, V> Forest<K, V> {
    pub fn new() -> Forest<K, V> {
        Forest {
            forest: HashMap::new(),
            orbits: HashSet::new()
        }
    }

    pub fn get_orbits<'a>(&'a self) -> impl Iterator<Item = &'a Set<V>> {
        self.orbits.iter()
            .filter_map(|key| {
                self.forest.get(key)
            })
            .filter_map(|e| {
                let ptr = unsafe { (&*e.ptr.as_ptr()).as_ref() };
                ptr.right()
            })
    }

    pub fn associate_with(&mut self, one: &K, two: &K, data: V) {
        match (self.forest.get(one), self.forest.get(two)) {
            (None, None) => {
                let op = SetPtr::new(Set::new(data));
                self.orbits.insert(one.clone());
                self.forest.insert(two.clone(), op.point());
                self.forest.insert(one.clone(), op);
            },
            (Some(x), None) => {
                self.forest.insert(two.clone(), x.point());
            },
            (None, Some(y)) => {
                self.forest.insert(one.clone(), y.point());
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

impl<K: Eq + Clone + std::hash::Hash, V: Default> Forest<K, V> {
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
        let mut forest: Forest<u32, ()> = Forest::new(); 
        let assocs = vec![(1, 2), (2, 3), (4, 5), (6, 7), (8, 9), (6, 2), (9, 4)];
        for (x,y) in assocs {
            forest.associate_with(&x, &y, ());
        }
        let orbits: Vec<&Set<()>> = forest.get_orbits().collect();
        assert_eq!(orbits.len(), 2);
    }
}
