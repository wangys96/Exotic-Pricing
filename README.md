Exotics Pricing
===
This repo has some applications of QuantLib.

### Projects
+ **AsianOption (C++)**   
 Implemented `DiscreteGeometricAverageStrikeEngine` and a pricing example
+ **Autocall CPU(C++&Pybind)**   
 Use `Pybind11` to wrap `PathGenerator` and `RandomSequenceGenerator` for fast Monte-Carlo simulation. An example of Autocall note is given.
 + **Autocall GPU(cuda&cupy)**   
 Use `cupy` and raw `cuda` kernel to perform Monte-Carlo simulation. An example of Autocall note is given and there is a 400x speedup compared with CPU.

