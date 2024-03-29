use libbgs_util::intpow;

/// An iterator yielding all of the factors of some number beneath a limit.
/// The type parameter `L` is the length of the factorization.
pub struct DivisorStream<'a> {
    source: &'a [(u128, usize)],
    stack: Vec<(usize, Vec<usize>)>,
    limit: u128,
    maximal_only: bool,
}

impl<'a> DivisorStream<'a> {
    /// Creates a new `DivisorStream`, which will return all of the factors of `source` beneath
    /// `limit`.
    /// In particular, it will return all values $d$ satisfying these properties:
    /// * $d | n$
    /// * $d < limit$
    /// * (if and only if `maximal_only` is True) There does not exist a $k$, $d | k | n$, with $k <
    /// limit$
    pub fn new(source: &'a [(u128, usize)], limit: u128, maximal_only: bool) -> DivisorStream {
        DivisorStream {
            source,
            limit,
            stack: vec![(0, vec![0; source.len()])],
            maximal_only,
        }
    }
}

impl<'a> Iterator for DivisorStream<'a> {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Vec<usize>> {
        let Some((i, state)) = self.stack.pop() else {
            return None;
        };
        // println!("{state:?}");
        let prod: u128 = state
            .iter()
            .zip(self.source)
            .map(|(d, (p, _))| intpow::<0>(*p, *d as u128))
            .product();
        let mut maximal = true;
        for j in i..self.source.len() {
            if state[j] == self.source[j].1 {
                continue;
            }
            // Assumption: the primes in a factorization are in increasing order.
            // If not, the break below should be a continue.
            if prod * self.source[j].0 > self.limit {
                break;
            }
            let mut next = state.clone();
            next[j] += 1;
            self.stack.push((j, next));
            maximal = false;
        }
        let Some((nonexhausted, _)) = self
            .source
            .iter()
            .enumerate()
            .find(|(j, (_, d))| state[*j] < *d)
        else {
            return Some(state);
        };
        maximal &= prod * self.source[nonexhausted].0 > self.limit;
        if self.maximal_only && !maximal {
            self.next()
        } else {
            Some(state)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stream_max() {
        let facts = [(2, 3), (3, 2), (5, 1)];
        let count = DivisorStream::new(&facts, 25, true).count();
        assert_eq!(count, 4);
    }

    #[test]
    fn test_stream_all() {
        let facts = [(2, 3), (3, 2), (5, 1)];
        let stream = DivisorStream::new(&facts, 25, false);
        let mut count = 0;
        for _ in stream {
            count += 1;
        }
        assert_eq!(count, 14);
    }

    #[test]
    fn test_stream_maximal_only() {
        let facts = [(2, 1), (7, 1), (13, 1)];
        let count = DivisorStream::new(&facts, 10, true).count();
        assert_eq!(count, 2);
    }
}
