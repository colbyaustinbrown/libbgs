`libbgs` is a collection of tools for investigating the structure of the Markoff graph modulo $p$.
* For the primary motivation and application of this library, see my preprint "An almost linear time algorithm testing whether the Markoff graph modulo $p$ is connected", available at [arxiv:2401.00630](https:://arxiv.org/abs/2401.00630).
* For an introduction to Markoff numbers, see Martin Aigner's "Markov's Theorem and 100 Years of the Uniqueness Conjecture".
* For an introduction to the Bourgain, Gamburd, and Sarnak algorithm, see their paper "Markoff Surfaces and Strong Approximation: 1", found at [arxiv:1607.01530](https://arxiv.org/abs/1607.01530).

Documentation in the form of rustdocs can be found on [my website](https://www.math.ucdavis.edu/~colbyabrown/libbgs-rustdocs/libbgs/index.html).

The examples directory contains a few sample programs using the `libbgs` algorithm, as well as the main entry point for the algorithm discussed in my preprint, [exhaustive-search.rs](examples/exhaustive-search.rs).
