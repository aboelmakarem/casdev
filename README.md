
CASDev Library

Ahmed M. Hussein

amhussein4@gmail.com

What is CASDev:
---------------
This repository holds the source code for the CASDev library. CASDev is 
a collection of tools that I have developed over the years and used as 
core components in solving various problems in academic and industrial 
research. The library components are developed as the need arises and do 
not follow any scheme or plan. Many of these components were written and 
re-written multiple times in order to simplify their design, extend their 
functionalities, optimize their performance or fix some of their bugs. 

CASDev was designed with two goals in mind, accuracy and performance. While 
more concise and more abstract code could have been written, this was 
sacrificed in order to develop something that is optimally geared towards 
correct execution, low run time and minimal memory footprint. 

The components of the library span a wide range of applications including 
numerical analysis, computational geometry, optimization algorithms, image 
and signal processing and machine learning models to name a few. Not all 
components are always available. I will add and remove components as I 
find bugs or work on their development. The latest commits represent the 
version I am most confident is bug free and optimal.New components will 
be added over time. Currently, I have a few components under development 
in the areas of machine learning algorithms, encryption and distributed 
computing. 

CASDev is actively developed, tested and maintained by Ahmed Hussein. The 
tests are not comprehensive so use it at your own risk. If you find a 
bug or have an idea on how to improve any of the components, please do 
not hesitate to reach me at amhussein4@gmail.com with your suggestions. 
I will do my best to address all comments and suggestions in a timely 
manner. 

How to build and use:
---------------------
CASDev is built as a shared object file that can be called from other 
applications. To build CASDev, run the python build script "build.py" 
which is located in the library's root directory. For an example on how 
to use the library in other applications, check out the build script of 
any of the test programs in the tests directory. 

Source Code:
------------
CASDev source code is distributed over four main directories. 

The "src" directory contains the implementation of all the library's 
functions. All functionality is written in C programming language. 

The "include" directory contains the header files to be included when 
linking against the library. 

The "tests" directory has a number of test programs that test the accuracy 
and performance of the different library components. All test code is 
written in C++ programming language. 

The "pymodules" directory has python modules that can be used to interface 
programs written in python to the library components. 

Why is CASDev implemented in C:
-------------------------------
My main development language is C++ and I use it almost every day in my 
work and research. I chose C for the library core components however for 
the following reasons

1. CASDev is designed for accuracy and performance, C is one of the fastest 
programming languages since it is very close to the machine compared to 
other high-level, memory-managed languages. 
2. All the core functionality does not need C++ specific features. None 
of the core modules are generic and there is no class heirarchy in place. 
3. Exporting C functionality is a lot easier than exporting C++ functionality 
since C compilers do not mangle names. Using the library with applications 
writtein in C++, python, or any other language is straightforward. 
4. C adds the least amount of overhead compared to any other programming 
language since there are no virtual table lookups, no function overloading 
or namespaces. 

This comes at the expense of the code being verbose and may be repeating 
accross components but all effort was made to reduce code redundancy. 


