//! Types and utilities for manipulating numbers in various types of finite fields.
//!
//! This module contains modular integers (i.e., $\mathbb{Z} / p\mathbb{Z}$ for prime $p$),
//! their quadratic finite field extensions (i.e., $\mathbb{Z} / p^2\mathbb{Z}$ for prime $p$), and
//! decompositions into direct sums of Sylow subgroups.
pub mod factor_stream;
pub mod factorization;
pub mod fp;
pub mod group;
pub mod quad_field;
pub mod sylow;
pub mod sylow_stream;
