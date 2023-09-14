/// A Bloom Filter, a probabilistic set.
/// Elements may be added to the filter, and then the filter may be tested for membership, with
/// false positives. The false positivity rate is determined by the size of the Bloom filter and
/// the number of hashes.
pub struct BloomFilter<T> {
    masks: Vec<u8>,
    hashes: Vec<Box<dyn Fn(&T) -> usize>>,
}

impl<T> BloomFilter<T> {
    /// Create a new Bloom filter, with the given size in bits and the given list of hashes to be
    /// applied to all members on addition and query.
    pub fn new(bits: usize, hashes: Vec<Box<dyn Fn(&T) -> usize>>) -> BloomFilter<T> {
        BloomFilter {
            masks: vec![0; bits >> 3],
            hashes,
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
    pub fn is_member<F>(&self, elem: &T, confirm: F) -> bool 
    where
        F: Fn(&T) -> bool,
    {
        self.is_member_prob(elem) && confirm(elem)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::*;

    #[test]
    fn test_bloom_filter() {
        let mut hashes = Vec::<Box<dyn Fn(&u128) -> usize>>::new();
        hashes.push(Box::new(|x| (x % 10_000) as usize));
        hashes.push(Box::new(|x| ((x >> 32) % 10_000) as usize));
        let mut filter = BloomFilter::<u128>::new(10_000, hashes);
        for i in 100_000..101_000 {
            filter.add(& intpow(i * 1000 + i * 10 + i, 2, 0));
        }
        for i in 100_000..100_500 {
            let x = intpow(i * 1000 + i * 10 + i, 2, 0);
            let check = filter.is_member_prob(&x);
            assert!(check);
        }
        let mut all = true;
        for i in 1_501..2_000 {
            let x = intpow(i * 1000 + i * 10 + i, 2, 0);
            all &= filter.is_member_prob(&x);
            if !all { break; } 
        }
        assert!(!all);
    }
}
