pub mod cli;
pub mod engine;
pub mod kmer;
pub mod peaks;
pub mod seqio;

pub use cli::{Config, parse_args};
pub use engine::{RunSummary, run};
