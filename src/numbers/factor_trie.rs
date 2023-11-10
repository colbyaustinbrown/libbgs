#[derive(Clone, Debug)]
pub struct FactorTrie<const L: usize, T> {
    i: usize,
    ds: [usize; L],
    children: [Option<Box<FactorTrie<L, T>>>; L],
    data: T,
}

#[derive(PartialEq, Eq, Clone, Copy)]
enum LeqMode {
    LEQ,
    STRICT,
}

pub trait FactorVisitor<const L: usize, T> {
    fn visit(&mut self, node: &FactorTrie<L, T>) {
        self.super_visit(node);
    }

    fn super_visit(&mut self, node: &FactorTrie<L, T>) {
        node.children.iter()
            .filter_map(|o| o.as_ref())
            .for_each(|n| {
                self.visit(n);
            });
    }
}

pub trait MutFactorVisitor<const L: usize, T> {
    fn visit(&mut self, node: &mut FactorTrie<L, T>) {
        self.super_visit(node);
    }

    fn super_visit(&mut self, node: &mut FactorTrie<L, T>) {
        node.children.each_mut()
            .into_iter()
            .filter_map(|o| o.as_mut())
            .for_each(|n| {
                self.visit(n);
            }); 
    }
}

impl<const L: usize, T> FactorTrie<L, T> {
    pub fn new(data: T) -> Self {
        FactorTrie {
            i: 0,
            ds: [0; L],
            data,
            children: std::array::from_fn(|_| None),
        }
    }

    pub fn get_or_new_child(&mut self, i: usize, data: T) -> &mut FactorTrie<L, T> {
        self.children[i].get_or_insert(Box::new(FactorTrie {
            i,
            ds: {
                let mut ds = self.ds;
                ds[i] += 1;
                ds
            },
            data,
            children: std::array::from_fn(|_| None),
        }))
    }

    pub fn add<F>(&mut self, t: [usize; L], gen: F)
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        self.add_helper(t, gen, LeqMode::STRICT);
    }

    pub fn add_leq<F>(&mut self, t: [usize; L], gen: F)
    where
        F: Fn(&[usize; L], usize) -> T,
    {
        self.add_helper(t, gen, LeqMode::LEQ);
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
                        children: std::array::from_fn(|_| None),
                    })
                })
                .add_helper(t, &gen, leq);
            if leq == LeqMode::STRICT { break; }
        }
    }

    pub fn map<S, F>(self, f: &mut F) -> FactorTrie<L, S>
    where
        F: FnMut(T, &[usize; L], usize) -> S,
    {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: f(self.data, &self.ds, self.i),
            children: self.children.map(|o| o.map(|n| Box::new(n.map(f)))),
        }
    }

    pub fn as_ref<'a>(&'a self) -> FactorTrie<L, &'a T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &self.data,
            children: self.children.each_ref().map(|o| {
                o.as_ref()
                    .map(|b| Box::new(FactorTrie::as_ref(Box::as_ref(b))))
            }),
        }
    }

    pub fn as_mut<'a>(&'a mut self) -> FactorTrie<L, &'a mut T> {
        FactorTrie {
            i: self.i,
            ds: self.ds,
            data: &mut self.data,
            children: self.children.each_mut().map(|o| {
                o.as_mut()
                    .map(|b| Box::new(FactorTrie::as_mut(Box::as_mut(b))))
            }),
        }
    }

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

    pub fn index(&self) -> usize {
        self.i
    }

    pub fn ds(&self) -> &[usize; L] {
        &self.ds 
    }

    pub fn data(&self) -> &T {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut T {
        &mut self.data
    }

    pub fn children(&self) -> &[Option<Box<FactorTrie<L, T>>>] {
        &self.children
    }
}
