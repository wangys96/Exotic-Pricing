/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <pybind11.h>
#include <Generator.h>

using namespace pybind11;


PYBIND11_MODULE(MCPath, m) {

    m.def("GeneratePath", &GeneratePath, "QuantLib QMC Path Generator",
          "today"_a, "num"_a, "steps"_a, "tenor"_a,
          "ir_type"_a,  "ir_term"_a,  "ir_data"_a,  "ir_dc"_a,
          "d_type"_a,   "d_term"_a,   "d_data"_a,   "d_dc"_a,
          "vol_type"_a, "vol_term"_a, "vol_data"_a, "vol_dc"_a,          
          "upout_type"_a,  "upout_ob"_a,   "upout_barrier"_a,
          "downout_type"_a,"downout_ob"_a, "downout_barrier"_a,
          "proc_type"_a, "output_matrix"_a,
          "bb"_a = true, "skip"_a = 0, "seed"_a = 42);

    m.def("GenerateRS", &GenerateRS, "QuantLib Sobol Random Seuqence Generator",
         "num"_a,  "steps"_a,  "tenor"_a,  "output_matrix"_a,  "bb"_a = true, "skip"_a = 0,"seed"_a = 42);

}