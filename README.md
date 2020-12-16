# pyflagsercount

The goal is to offer some python binding to flagser c++ code, with some added features.
- compute simplex count over a filtration
- compute maximal simplex count
- return list of simplex as a simplex tree(one tree for each vertices)

Project derived from Jason P. Smith pyflagsercontain (itself derived from Daniel Luetgehetmann flagser (itself derived from ripser))


## Installation from source

You will need  
- cmake 
- google sparse hash
- pybind11 as a subdirectory

Do at the root of the project:  
`$ git clone https://github.com/sparsehash/sparsehash`  
You have to then compile it, no need to install:  
`$ cd sparsehash && ./configure && make && cd ..  `

`$ git clone https://github.com/pybind/pybind11.git `

pybind11 is header only so no compilation is needed  

## Installation from wheel

In dist/ you have some wheel package that may be suitable for your system. You can try to install it this way:  
`$ pip install pyflagsercount-0.0.1-cp36-cp36m-linux_x86_64.whl `  
For linux for example 