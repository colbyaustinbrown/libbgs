use crate::util::intpow;

/// An iterator yielding all of the factors of some number beneath a limit.
/// The type parameter `L` is the length of the factorization.
pub struct FactorStream<'a, const L: usize> {
    source: &'a [(u128, usize)],
    stack: Vec<(usize, [usize; L])>,
    limit: u128,
    maximal_only: bool,
}

impl<'a, const L: usize> FactorStream<'a, L> {
    /// Creates a new `FactorStream`, which will return all of the factors of `source` beneath
    /// `limit`.
    /// In particular, it will return all values $d$ satisfying these properties:
    /// * $d | n$
    /// * $d < limit$
    /// * (if and only if `maximal_only` is True) There does not exist a $k$, $d | k | n$, with $k <
    /// limit$
    pub fn new(source: &'a [(u128, usize)], limit: u128, maximal_only: bool) -> FactorStream<L> {
        FactorStream {
            source,
            limit,
            stack: vec![(0, [0; L])],
            maximal_only,
        }
    }
}

impl<'a, const L: usize> Iterator for FactorStream<'a, L> {
    type Item = [usize; L];

    fn next(&mut self) -> Option<[usize; L]> {
        let Some((i, state)) = self.stack.pop() else { return None; };
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
            let mut next = state;
            next[j] += 1;
            self.stack.push((j, next));
            maximal = false;
        }
        maximal &= prod * self.source[0].0 > self.limit;
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
        let stream = FactorStream::<3>::new(&facts, 25, true);
        let mut count = 0;
        for _ in stream {
            count += 1;
        }
        assert_eq!(count, 4);
    }

    #[test]
    fn test_stream_all() {
        let facts = [(2, 3), (3, 2), (5, 1)];
        let stream = FactorStream::<3>::new(&facts, 25, false);
        let mut count = 0;
        for _ in stream {
            count += 1;
        }
        assert_eq!(count, 14);
    }
}
