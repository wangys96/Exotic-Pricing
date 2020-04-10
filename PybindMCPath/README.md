PyBind Monte Carlo Path Generator
===
This project uses pybind11 to implement the function generating MC paths through QuantLib.   

### Compile  
To complile `MCPath`, remember to set `IncludePath` with `.../pybind11`, `.../<PythonPath>/include`, `.../<QuantLibPath>/ql`; and set `LibraryPath` with `.../<PythonPath>/libs`, `.../<QuantLibPath>/libs`.  

### Performance
```shell
PS C:\Apps\QuantLib-1.18\Examples\MCPath\bin> python test.py
Test generating MC paths...
 [Result]:  27.41994833946228
Test generating MC paths with early stop barrier ...
 [Result]:  10.283002376556396
Test generating MC paths with 4 processes...
 [Result]:  12.773000478744507
Test generating MC paths with 4 processes and early stop barrier...
 [Result]:  6.103964328765869
```

### Usage
```python
import MCPath
MCPath.GeneratePath(
    today: tuple,                         # (day,month,year)
    num: int,                             # number of simulation paths
    steps: int,                           # number of steps (e.g. 366 days)
    tenor: float,                         # time length in unit of year (e.g. 366/365 year) 
    
    ir_type: int,                         # interest rate type, 0=flat, 1=spot, 2=forward, 3=discountfactor
    ir_term: numpy.ndarrayint32,          # days after today
    ir_data: numpy.ndarrayfloat64,        # rate at each term
    ir_dc: int,                           # daycounter, 0=Act/365, 1=Act/Act, 2=Act/360, 3=30/360
    
    d_type: int,                          # dividend rate type, same as ir_type
    d_term: numpy.ndarrayint32,           # days after today
    d_data: numpy.ndarrayfloat64,         # rate at each term
    d_dc: int,                            # daycounter, same as ir_dc
    
    vol_type: int,                        # volatility type, 0=flat, 1=spot
    vol_term: numpy.ndarrayint32,         # days after today
    vol_data: numpy.ndarrayfloat64,       # vol rate at each term
    vol_dc: int,                          # daycounter, same as ir_dc
    
    upout_type: int,                      # up knock-out type, 0=NoBarrier, 1=ConstBarrier, 2=NonConstBarrier
    upout_ob: numpy.ndarraybool,          # boolean array that if is knock observation days, length=steps+1
    upout_barrier: numpy.ndarrayfloat64,  # barrier value, if type=1 only 1 number is needed in this numpy array

    downout_type: int,                    # down knock-out type, same as upout_type
    downout_ob: numpy.ndarraybool,        # boolean array, same as upout_ob
    downout_barrier: numpy.ndarrayfloat64,# barrier value, same as upout_barrier
    proc_type: int,                       # type of stochastic process, 1=BS, 2=BSM(with dividend)
    input_matrix: numpy.ndarrayfloat64,   # an empty numpy array with shape(num,steps+1)
    bb: bool = True,                      # use Brownian Bridge
    skip: int = 0,                        # skip random sequence
    seed: int = 42
)
```
