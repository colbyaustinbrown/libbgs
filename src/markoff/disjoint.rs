use std::collections::{HashSet, HashMap};
use std::rc::Rc;
use std::cmp::Ordering;
use either::{Left, Right, Either};
use std::cell::{Ref, RefCell};

#[derive(Debug)]
pub struct Set<V> {
    rank: u16,
    pub data: V
}

#[derive(Debug)]
pub struct SetPtr<V> {
    ptr: Rc<RefCell<Either<SetPtr<V>, Set<V>>>>
}

// TODO: disjoint should not be pub (nor should SetPtr above)
// made public for debugging purposes
pub struct Disjoint<K, V> {
    default: V,
    combine: Box<dyn Fn(&V, &V) -> V>,
    pub disjoint: HashMap<K, SetPtr<V>>,
    orbits: HashSet<K>
}

impl<K, V> Disjoint<K, V> 
where
    K: Eq + Clone + std::hash::Hash,
    V: Clone 
{
    pub fn new<F>(default: V, combine: F) -> Disjoint<K,V>
    where
        for<'a> F: Fn(&V, &V) -> V + 'a
    {
        Disjoint {
            default,
            combine: Box::new(combine),
            disjoint: HashMap::new(),
            orbits: HashSet::new()
        }
    }

    pub fn get_orbits(&self) -> impl Iterator<Item = (&K, Ref<Set<V>>)> {
        self.orbits.iter()
            .filter_map(|key| {
                self.disjoint.get(key).map(|r| (key, r))
            })
            .filter_map(|(key, op)| {
                Ref::filter_map(op.ptr.borrow(), |e| {
                    e.as_ref().right()
                }).ok().map(|r| (key, r))
            })
    }

    pub fn associate(&mut self, one: &K, two: &K) {
        match (self.disjoint.get(one), self.disjoint.get(two)) {
            (None, None) => {
                let op = SetPtr::new(Set::new(self.default.clone()));
                self.orbits.insert(one.clone());
                self.disjoint.insert(two.clone(), op.point());
                self.disjoint.insert(one.clone(), op);
            },
            (Some(x), None) => {
                let (root,_) = unsafe { x.root() };
                self.disjoint.insert(two.clone(), root.point());
            },
            (None, Some(y)) => {
                let (root,_) = unsafe { y.root() };
                self.disjoint.insert(one.clone(), root.point());
            },
            (Some(x), Some(y)) => {
                let (p1,o1) = unsafe { x.root() };
                let (p2,o2) = unsafe { y.root() };
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
                            data: (self.combine)(&o1.data, &o2.data)
                        }));
                    }
                }
            }
        }
    }

    pub fn update(&mut self, k: &K, v: V) {
        let Some(orbit) = self.disjoint.get(k) else { return; };
        let (ptr,set) = unsafe { orbit.root() };
        ptr.ptr.replace(Right(Set {
            rank: set.rank,
            data: v
        }));
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

    unsafe fn root(&self) -> (&SetPtr<V>, &Set<V>) {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assoc() {
        let mut disjoint: Disjoint<u32, ()> = Disjoint::new((), |_,_| ()); 
        let assocs = vec![(1, 2), (2, 3), (4, 5), (6, 7), (8, 9), (6, 2), (9, 4)];
        for (x,y) in assocs {
            disjoint.associate(&x, &y);
        }
        let orbits: Vec<(&u32, Ref<Set<()>>)> = disjoint.get_orbits().collect();
        assert_eq!(orbits.len(), 2);
    }
}
