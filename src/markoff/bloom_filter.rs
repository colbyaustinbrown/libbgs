pub struct BloomFilter<T> {
    masks: Vec<bool>,
    hashes: Vec<Box<dyn Fn(&T) -> usize>>,
}

impl<T> BloomFilter<T> {
    pub fn new(bits: usize, hashes: Vec<Box<dyn Fn(&T) -> usize>>) -> BloomFilter<T> {
        BloomFilter {
            masks: vec![false; bits],
            hashes,
        }
    }

    pub fn add(&mut self, elem: &T) {
        self.hashes.iter()
            .for_each(|hash| {
                self.masks[hash(elem)] = true;
            });
    }

    pub fn is_member_prob(&self, elem: &T) -> bool {
        self.hashes.iter()
            .all(|hash| {
                self.masks[hash(elem)]
            })
    }

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
