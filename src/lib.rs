// #![feature(array_methods)]
#![warn(missing_docs, unused_lifetimes, unused_tuple_struct_fields, variant_size_differences, clippy::all)]
#![allow(clippy::new_without_default, clippy::needless_range_loop)]

//! A library for the manipulation and investigation of Markoff numbers and Markoff graphs modulo
//! primes, specifically with a focus on the Bourgain, Gamburd, and Sarnak algorithm.
pub mod markoff;
pub mod numbers;
pub mod streams;
