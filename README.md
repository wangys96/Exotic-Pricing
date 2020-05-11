Exotics Pricing
===
This repo store some applications of QuantLib.

### Projects
+ **AsianOption (C++)**   
 Implemented `DiscreteGeometricAverageStrikeEngine` and a pricing example
+ **Pybind+QuantLib Montecarlo (C++&Python)**   
 Use `Pybind11` to wrap `PathGenerator` and `RandomSequenceGenerator` for fast Monte-Carlo simulation. An example of Autocall note is given.
 + **Cuda+Cupy MonteCarlo (C++&Python)**   
 Use `cupy` and raw `cuda` kernel to perform Monte-Carlo simulation. An example of Autocall note is given and there is a 400x speedup compared with CPU.

