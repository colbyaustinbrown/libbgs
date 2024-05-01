#![warn(missing_docs, unused_lifetimes, dead_code, variant_size_differences, clippy::all)]
#![allow(clippy::new_without_default, clippy::needless_range_loop)]

//! A library for the manipulation and investigation of Markoff numbers and Markoff graphs modulo
//! primes, specifically with a focus on the Bourgain, Gamburd, and Sarnak algorithm.
pub mod markoff;
pub mod numbers;
pub mod streams;
pub mod util;
