use crate::util::intpow;

pub struct FactorStream<'a> {
    source: &'a [(u128, u128)],
    stack: Vec<(usize, Vec<u128>)>,
    limit: u128,
    maximal_only: bool,
}

impl<'a> FactorStream<'a> {
    pub fn new(source: &'a [(u128, u128)], limit: u128) -> FactorStream {
        FactorStream {
            source,
            limit,
            stack: vec![(0, vec![0; source.len()])],
            maximal_only: true,
        }
    }
}

impl<'a> Iterator for FactorStream<'a> {
    type Item = Vec<u128>;

    fn next(&mut self) -> Option<Vec<u128>> {
        let Some((i, state)) = self.stack.pop() else { return None; };
        // println!("{state:?}");
        let prod: u128 = state
            .iter()
            .zip(self.source)
            .map(|(d, (p, _))| intpow(*p, *d, 0))
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
    fn test_stream() {
        let facts = vec![(2, 3), (3, 2), (5, 1)];
        let stream = FactorStream::new(&facts, 25);
        let mut count = 0;
        for x in stream {
            println!("{x:?}");
            // assert!(x < 25);
            count += 1;
        }
        assert_eq!(count, 4);
    }
}
