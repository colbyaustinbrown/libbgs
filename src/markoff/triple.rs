use crate::numbers::FpNum;

/// A Markoff triple modulo `P`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Triple<const P: u128>(FpNum<P>, FpNum<P>, FpNum<P>);

/// The three coordinates of a Markoff triple.
pub enum Pos {
    /// The first coordinate.
    A,
    /// The second coordinate.
    B,
    /// The third coordinate.
    C,
}

impl<const P: u128> Triple<P> {
    /// Creates a Markoff triple, but only if `a`, `b`, and `c` are a solution.
    pub fn new(a: FpNum<P>, b: FpNum<P>, c: FpNum<P>) -> Option<Triple<P>> {
        if a * a + b * b + c * c - a * b * c == FpNum::from(0) {
            Some(Triple(a, b, c))
        } else {
            None
        }
    }

    /// Creates a new Markoff triple, without checking that it is actually a solution.
    pub fn new_unchecked(a: FpNum<P>, b: FpNum<P>, c: FpNum<P>) -> Triple<P> {
        Triple(a, b, c)
    }

    /// Returns the Markoff triple obtained by applying the Vieta involution to the Markoff triple.
    pub fn vieta(&self, coord: Pos) -> Triple<P> {
        match coord {
            Pos::A => Triple(self.1 * self.2 - self.0, self.1, self.2),
            Pos::B => Triple(self.0, self.0 * self.2 - self.1, self.2),
            Pos::C => Triple(self.0, self.1, self.0 * self.1 - self.2),
        }
    }

    /// Returns the first coordinate.
    #[inline(always)]
    pub fn a(&self) -> FpNum<P> {
        self.0
    }

    /// Returns the second coordinate.
    #[inline(always)]
    pub fn b(&self) -> FpNum<P> {
        self.1
    }

    /// Returns the third coordinate.
    #[inline(always)]
    pub fn c(&self) -> FpNum<P> {
        self.2
    }

    /// Returns the designated coordinate.
    #[inline(always)]
    pub fn get(&self, p: Pos) -> FpNum<P> {
        match p {
            Pos::A => self.0,
            Pos::B => self.1,
            Pos::C => self.2,
        }
    }
}
