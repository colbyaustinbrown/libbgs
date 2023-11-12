use crate::numbers::Factorization;
use std::sync::Arc;

#[derive(Clone, Debug)]
/// A trie of prime factors in increasing order; that is, a none with word $p$ will have
/// only children with word $q \geq p$.
pub struct FactorTrie<const L: usize, T> {
    i: usize,
    ds: [usize; L],
    fact: Arc<Factorization>,
    children: [Option<Box<FactorTrie<L, T>>>; L],
    /// Data associated with the key given by the concatenation of this node's ancestors' words.
    pub data: T,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LeqMode {
    Leq,
    Strict,
}

impl<const L: usize, T> FactorTrie<L, T> {
    /// Creates a new trie associated to the given `factorization`.
    /// The trie begins as only a root node containing `data`; children must be explicitly added
    /// via the `[add]`, `[add_leq]`, or `[get_or_new_child]` methods.
    pub fn new(factorization: Factorization, data: T) -> Self {
        FactorTrie {
            i: 0,
            ds: [0; L],
            data,
            fact: Arc::new(factorization),
            children: std::array::from_fn(|_| None),
        }
    }

    /// Returns this node's child at index `i`, or creates the child, initialized with the result
    /// of the lazily-evaluated `data`.
    pub fn get_or_new_child<F>(&mut self, i: usize, data: F) -> &mut FactorTrie<L, T> 
    where
        F: FnOnce() -> T,
    {
        self.children[i].get_or_insert(Box::new(FactorTrie {
            i,
            ds: {
                let mut ds = self.ds;
                ds[i] += 1;
                ds
            },
            data: data(),
            fact: Arc::clone(&self.fact),
            children: std::array::from_fn(|_| None),
        }))
    }

    /// Add a new child to the position representated by `t`, creating all necessary new children
    /// along the way via `gen`.
    pub fn add<F>(&mut self, t: [usize; L], gen: F)
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        self.add_helper(t, gen, LeqMode::Strict);
    }

    /// Add a new child and all children represented by `t` and all its divisors (subordinate
    /// arrays), creating all necessary new children along the way via `gen`.
    pub fn add_leq<F>(&mut self, t: [usize; L], gen: F)
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        self.add_helper(t, gen, LeqMode::Leq);
    }

    fn add_helper<F>(&mut self, t: [usize; L], gen: F, leq: LeqMode)
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        for j in self.i..L {
            if t[j] <= self.ds[j] {
                continue;
            }
            self.children[j]
                .get_or_insert_with(|| {
                    let mut ds = self.ds;
                    ds[j] += 1;
                    Box::new(FactorTrie {
                        i: j,
                        ds,
                        data: gen(&ds, j),
                        fact: Arc::clone(&self.fact),
                        children: std::array::from_fn(|_| None),
                    })
                })
                .add_helper(t, &gen, leq);
            if leq == LeqMode::Strict { break; }
        }
    }

    /// Transforms this trie into an equivalent trie with the same shape, but all data mapped via
    /// `f`.
    pub fn map<S, F>(self, f: &mut F) -> FactorTrie<L, S>
    where
        F: FnMut(T, &[usize; L], usize) -> S,
    {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: f(self.data, &self.ds, self.i),
            fact: Arc::clone(&self.fact),
            children: self.children.map(|o| o.map(|n| Box::new(n.map(f)))),
        }
    }

    /// Returns a trie of borrowed data.
    pub fn as_ref(&self) -> FactorTrie<L, &T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &self.data,
            fact: Arc::clone(&self.fact),
            children: self.children.each_ref().map(|o| {
                o.as_ref()
                    .map(|b| Box::new(FactorTrie::as_ref(Box::as_ref(b))))
            }),
        }
    }

    /// Returns a trie of mutably borrowed data.
    pub fn as_mut(&mut self) -> FactorTrie<L, &mut T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &mut self.data,
            fact: Arc::clone(&self.fact),
            children: self.children.each_mut().map(|o| {
                o.as_mut()
                    .map(|b| Box::new(FactorTrie::as_mut(Box::as_mut(b))))
            }),
        }
    }

    /// Runs `f` on each node, in a pre-order traversal.
    pub fn for_each<F>(&self, f: &mut F)
    where
        F: FnMut(&T, [usize; L])
    {
        f(&self.data, self.ds);
        self.children.iter()
            .for_each(|o| match o {
                None => {},
                Some(b) => b.for_each(f),
            });
    }

    /// The index of this node's word in the prime factorization array.
    pub fn index(&self) -> usize {
        self.i
    }

    /// The array of powers on this node's prime factoraziton.
    pub fn ds(&self) -> &[usize; L] {
        &self.ds 
    }

    /// The prime factorization represented by this trie.
    pub fn fact(&self) -> &Factorization {
        &self.fact
    }

    /// This node's array of children.
    pub fn children(&self) -> &[Option<Box<FactorTrie<L, T>>>] {
        &self.children
    }

    /// Returns a reference to the child at index `i`, if there is one.
    pub fn child(&self, i: usize) -> Option<&FactorTrie<L, T>> {
        self.children[i].as_deref()
    }

    /// Returns a mutable reference to the child at index `i`, if there is one.
    pub fn child_mut(&mut self, i: usize) -> Option<&mut FactorTrie<L, T>> {
        self.children[i].as_deref_mut()
    }
}
