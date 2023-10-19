//! Types and utilities for manipulating numbers in various types of finite fields.
//!
//! This module contains modular integers (i.e., $\mathbb{Z} / p\mathbb{Z}$ for prime $p$),
//! their quadratic finite field extensions (i.e., $\mathbb{Z} / p^2\mathbb{Z}$ for prime $p$), and
//! decompositions into direct sums of Sylow subgroups.
mod factorization;
mod fp;
mod group;
mod quad_field;
mod sylow;
mod montgomery;

pub use factorization::*;
pub use fp::*;
pub use group::*;
pub use quad_field::*;
pub use sylow::*;
pub use montgomery::*;
