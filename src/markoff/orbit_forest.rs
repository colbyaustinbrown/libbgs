use std::collections::{HashSet, HashMap};
use std::rc::Rc;
use std::cmp::Ordering;
use either::{Left, Right, Either};
use std::cell::{RefCell};

#[derive(Debug)]
struct Orbit {
    rank: u16,
    has_large: bool
}

#[derive(Debug)]
struct OrbitPtr {
    ptr: Rc<RefCell<Either<OrbitPtr, Orbit>>>
}

#[derive(Debug)]
struct Forest<K> {
    forest: HashMap<K, OrbitPtr>,
    orbits: HashSet<K>
}

impl<K: Eq + Clone + std::hash::Hash> Forest<K> {
    pub fn new() -> Forest<K> {
        Forest {
            forest: HashMap::new(),
            orbits: HashSet::new()
        }
    }

    pub fn get_orbits<'a>(&'a self) -> impl Iterator<Item = &'a Orbit> {
        self.orbits.iter()
            .filter_map(|key| {
                self.forest.get(key)
            })
            .filter_map(|e| {
                let ptr = unsafe { (&*e.ptr.as_ptr()).as_ref() };
                ptr.right()
            })
    }

    pub fn associate(&mut self, one: &K, two: &K) {
        match (self.forest.get(one), self.forest.get(two)) {
            (None, None) => {
                let op = OrbitPtr::new(Orbit::new());
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
                        p2.ptr.replace(Right(Orbit {
                            rank: o2.rank,
                            has_large: o2.has_large
                        }));
                    }
                }
            }
        }
    }
}

impl Orbit {
    fn new() -> Orbit {
        Orbit {
            rank: 0,
            has_large: false
        }
    }
}

impl OrbitPtr {
    fn new(orbit: Orbit) -> OrbitPtr {
        OrbitPtr { 
            ptr: Rc::new(RefCell::new(Right(orbit)))
        }
    }

    fn point(&self) -> OrbitPtr {
        OrbitPtr {
            ptr: Rc::clone(&self.ptr)
        }
    }

    fn root(&self) -> (&OrbitPtr, &Orbit) {
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
        let mut forest: Forest<u32> = Forest::new(); 
        let assocs = vec![(1, 2), (2, 3), (4, 5), (6, 7), (8, 9), (6, 2), (9, 4)];
        for (x,y) in assocs {
            forest.associate(&x, &y);
        }
        let orbits: Vec<&Orbit> = forest.get_orbits().collect();
        assert_eq!(orbits.len(), 2);
    }
}
