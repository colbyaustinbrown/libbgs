use std::marker::PhantomData;
use std::sync::Arc;

/// A Bloom Filter, a probabilistic set.
/// Elements may be added to the filter, and then the filter may be tested for membership, with
/// false positives. The false positivity rate is determined by the size of the Bloom filter and
/// the number of hashes.
pub struct BloomFilter<T, F> {
    masks: Vec<u8>,
    hashes: Arc<Vec<F>>,
    _phantom: PhantomData<T>,
}

impl<T, F> BloomFilter<T, F> 
where
    F: Fn(&T) -> usize + Send + Sync,
{
    /// Create a new Bloom filter, with the given size in bits and the given list of hashes to be
    /// applied to all members on addition and query.
    pub fn new(bits: usize, hashes: Vec<F>) -> BloomFilter<T, F> 
    {
        BloomFilter {
            masks: vec![0; bits >> 3],
            hashes: Arc::new(hashes),
            _phantom: PhantomData,
        }
    }

    /// Add `elem` to the Bloom filter.
    pub fn add(&mut self, elem: &T) {
        self.hashes.iter()
            .for_each(|hash| {
                let h = hash(elem);
                self.masks[h >> 3] |= 1 << (h & 0b111);
            });
    }

    /// True if `elem` is in the set.
    /// If `elem` is not in the set, this method returns False; i.e., this method return false
    /// positives, but not false negatives.
    pub fn is_member_prob(&self, elem: &T) -> bool {
        self.hashes.iter()
            .all(|hash| {
                let h = hash(elem);
                self.masks[h >> 3] & (1 << (h & 0b111)) != 0
            })
    }

    /// True if `elem` is in the set, lazily confirming the result with the `confirm` closure to
    /// guard against false positives.
    pub fn is_member<G>(&self, elem: &T, confirm: G) -> bool 
    where
        G: Fn(&T) -> bool,
    {
        self.is_member_prob(elem) && confirm(elem)
    }

    /// Modifies `self` to include elements from `other`.
    /// The false positivity rate of the resultant bloom filter will be greater than or equal to
    /// the maximum of the false positivity rates of the two operands.
    pub fn union(&mut self, other: &Self) {
        let l = usize::max(self.masks.len(), other.masks.len());
        for i in 0..l {
            self.masks[i] |= other.masks[i];
        }
    }
}

impl<T, F> Clone for BloomFilter<T, F> {
    fn clone(&self) -> BloomFilter<T, F> {
        BloomFilter {
            masks: self.masks.clone(),
            hashes: Arc::clone(&self.hashes),
            _phantom: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::*;

    #[test]
    fn test_bloom_filter() {
        let mut hashes = Vec::<Box<dyn Fn(&u128) -> usize + Send + Sync>>::new();
        hashes.push(Box::new(|x| (x % 10_000) as usize));
        hashes.push(Box::new(|x| ((x >> 32) % 10_000) as usize));
        let mut filter = BloomFilter::<u128, _>::new(10_000, hashes);
        for i in 100_000..101_000 {
            filter.add(& intpow::<0>(i * 1000 + i * 10 + i, 2));
        }
        for i in 100_000..100_500 {
            let x = intpow::<0>(i * 1000 + i * 10 + i, 2);
            let check = filter.is_member_prob(&x);
            assert!(check);
        }
        let mut all = true;
        for i in 1_501..2_000 {
            let x = intpow::<0>(i * 1000 + i * 10 + i, 2);
            all &= filter.is_member_prob(&x);
            if !all { break; } 
        }
        assert!(!all);
    }
}
