use std::marker::PhantomData;

use crate::numbers::FactoredSize;

type Child<S, const L: usize, C, T> = Box<FactorTrie<S, L, C, T>>;

/// A trie of prime factors in increasing order; that is, a none with word $p$ will have
/// only children with word $q \geq p$.
pub struct FactorTrie<S, const L: usize, C, T> {
    i: usize,
    ds: [usize; L],
    children: [Option<Child<S, L, C, T>>; L],
    /// Data associated with the key given by the concatenation of this node's ancestors' words.
    pub data: T,
    _phantom: PhantomData<(S, C)>,
}

impl<S, const L: usize, C> FactorTrie<S, L, C, ()>
where
    C: FactoredSize<S>,
{
    /// Creates a new trie associated to the given `factorization`.
    /// The trie contains the divisors of `C::FACTORS`.
    /// The words of a divisor are its prime factors, with multiplicities, in ascending order.
    pub fn new() -> FactorTrie<S, L, C, ()> {
        FactorTrie::new_with(|_, _| ())
    }
}

impl<S, const L: usize, C, T> FactorTrie<S, L, C, T> 
where
    C: FactoredSize<S>,
{
    /// Creates a new trie associated to the given `factorization`.
    /// The trie contains the divisors of `C::FACTORS`.
    /// The words of a divisor are its prime factors, with multiplicities, in ascending order.
    /// The trie is seeded with data generated by `f`.
    pub fn new_with<F>(f: F) -> FactorTrie<S, L, C, T> 
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        let mut res = FactorTrie {
            i: 0,
            ds: [0; L],
            data: f(&[0; L], 0),
            children: std::array::from_fn(|_| None),
            _phantom: PhantomData,
        };
        res.new_helper(std::array::from_fn(|i| C::FACTORS[i].1), &f);
        res
    }

    fn new_helper<F>(&mut self, t: [usize; L], f: &F) 
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        for j in self.i..L {
            if self.ds[j] >= t[j] {
                continue;
            }
            self.children[j]
                .get_or_insert_with(|| {
                    let mut ds = self.ds;
                    ds[j] += 1;
                    let mut child = FactorTrie {
                        i: j,
                        ds,
                        data: f(&ds, j),
                        children: std::array::from_fn(|_| None),
                        _phantom: PhantomData,
                    };
                    child.new_helper(t, f);
                    Box::new(child)
                });
        }
    }
}

impl<S, const L: usize, C, T> FactorTrie<S, L, C, T> {
    /// Returns this node's child at index `i`, or creates the child, initialized with the result
    /// of the lazily-evaluated `data`.
    pub fn get_or_new_child<F>(&mut self, i: usize, data: F) -> &mut FactorTrie<S, L, C, T>
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
            children: std::array::from_fn(|_| None),
            _phantom: PhantomData,
        }))
    }

    /// Updates a single node in the trie.
    pub fn update<F>(&mut self, t: &[usize; L], gen: F)
    where
        F: Fn(&[usize; L], &mut T)
    {
        for j in self.i..L {
            if t[self.i] == self.ds[self.i] { continue; }
            if let Some(ref mut child) = self.children[j] {
                child.update(t, gen);
                return;
            }
        }
        gen(&self.ds, &mut self.data);
    }

    /// Transforms this trie into an equivalent trie with the same shape, but all data mapped via
    /// `f`.
    pub fn map<U, F>(self, f: &F) -> FactorTrie<S, L, C, U>
    where
        F: Fn(T, &[usize; L], usize) -> U,
    {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: f(self.data, &self.ds, self.i),
            children: self.children.map(|o| o.map(|n| Box::new(n.map(f)))),
            _phantom: PhantomData,
        }
    }

    /// Returns a trie of borrowed data.
    pub fn as_ref(&self) -> FactorTrie<S, L, C, &T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &self.data,
            children: self.children.each_ref().map(|o| {
                o.as_ref()
                    .map(|b| Box::new(FactorTrie::as_ref(Box::as_ref(b))))
            }),
            _phantom: PhantomData,
        }
    }

    /// Returns a trie of mutably borrowed data.
    pub fn as_mut(&mut self) -> FactorTrie<S, L, C, &mut T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &mut self.data,
            children: self.children.each_mut().map(|o| {
                o.as_mut()
                    .map(|b| Box::new(FactorTrie::as_mut(Box::as_mut(b))))
            }),
            _phantom: PhantomData,
        }
    }

    /// Runs `f` on each node, in a pre-order traversal.
    pub fn for_each<F>(&self, f: &mut F)
    where
        F: FnMut(&T, [usize; L]),
    {
        f(&self.data, self.ds);
        self.children.iter().for_each(|o| match o {
            None => {}
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

    /// This node's array of children.
    pub fn children(&self) -> &[Option<Child<S, L, C, T>>] {
        &self.children
    }

    /// Returns a reference to the child at index `i`, if there is one.
    pub fn child(&self, i: usize) -> Option<&FactorTrie<S, L, C, T>> {
        self.children[i].as_deref()
    }

    /// Returns a mutable reference to the child at index `i`, if there is one.
    pub fn child_mut(&mut self, i: usize) -> Option<&mut FactorTrie<S, L, C, T>> {
        self.children[i].as_deref_mut()
    }
}

impl<S, const L: usize, C, T: Clone> Clone for FactorTrie<S, L, C, T> {
    fn clone(&self) -> Self {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: self.data.clone(),
            children: self.children.clone(),
            _phantom: PhantomData,
        }
    }
}
