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
git clone https://github.com/sparsehash/sparsehash  
compile it, no need to install: 
cd sparsehash && ./configure && make && cd ..  

git clone https://github.com/pybind/pybind11.git  

pybind11 is header only so no compilation is needed  