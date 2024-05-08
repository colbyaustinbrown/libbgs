//! Types and utilities for manipulating numbers in various types of finite fields.
//!
//! This module contains modular integers (i.e., $\mathbb{Z} / p\mathbb{Z}$ for prime $p$),
//! their quadratic finite field extensions (i.e., $\mathbb{Z} / p^2\mathbb{Z}$ for prime $p$), and
//! decompositions into direct sums of Sylow subgroups.
extern crate libbgs_macros;
mod factor_trie;
mod factorization;
mod fp;
mod group;
mod quad_num;
mod sylow;
mod norm1;

pub use factor_trie::*;
pub use factorization::*;
pub use fp::*;
pub use group::*;
pub use libbgs_macros::*;
pub use quad_num::*;
pub use sylow::*;
pub use norm1::*;
