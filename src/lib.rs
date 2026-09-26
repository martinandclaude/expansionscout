//! ExpansionScout: targeted repeat-expansion calling from ONT BAMs, with
//! methylation.
//!
//! A port of the Python reference implementation, written for
//! reproducibility rather than speed: one statically linkable binary, with
//! the STRchive release, the clinical overlay and the VCF dialects compiled
//! in, so that a binary and its version are enough to reproduce a result.
//!
//! The Python is developed in its own repository,
//! `martinandclaude/expansionscout-dev`, and remains the reference: its
//! parity tests run both on the same inputs and require identical output,
//! against the commit of this repository it pins. Paths in comments here
//! such as `expansionscout/engine.py` or `tests/test_engine.py` are files
//! there, and "the Python" is that implementation.

pub mod alleles;
pub mod bam;
pub mod catalog;
pub mod classify;
pub mod cli;
pub mod consensus;
pub mod decompose;
pub mod engine;
pub mod methylation;
pub mod npstat;
pub mod pyfmt;
pub mod qc;
pub mod reads;
pub mod report;
pub mod vcf;
