pub struct Factorization {
    pub value: u128,
    pub factors: Vec<u128>,
    pub primepowers: Vec<(u128, u128)>,
}

impl Factorization {
    pub fn len(&self) -> usize {
        self.primepowers.len()
    }
}

