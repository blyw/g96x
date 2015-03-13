# Project Description #

This project is all about _**exploring**_ ways to process molecular dynamic simulation trajectories in g96 format. Additionally, programs able to process the resulting data, either visualization or reanalysis, can be included.

## Developer guidelines ##
  * Each program has its own directory. It should be possible to compile each program independently.
  * All programs share a similar XML input parameter file construct.
  * All programs read and write data in compressed format (GZip). Uncompressed data is allowed for files with file size smaller than 20 MB.
  * Regardless of platform dependent performance differences, all programs can be compiled on Windows and major Linux distributions.
  * C++11 complaint code and C# is preferred over C.
  * C++11 complaint code is preferred over C#.
  * Python, JAVA and/or FORTRAN can be used for proof-of-concept, provided that the algorithm will be made  available in C++ or C# within reasonable time after proof-of-concept.
  * Python scripts that utilize matplotlib might be provided as extras and should not be considered as an integral part of any program.
  * Minimal dependencies on third party libraries. Especially those that are not actively maintained or developed. In case any of such dependency is abandoned by its community and developers, migration to an alternative should be done within reasonable time.
  * Extensive documentation and comments on algorithms. Inexperienced programmers should not have to puzzle to understand the code and/or the algorithm nor run into pitfalls that can only be identified by experienced users. Documentation and comments should be as easy as reading a cook book.
  * Clear error messages. 'Segmentation fault' is not a clear error message.
  * Programs should be suitable for inexperienced users, either due to lack of extensive knowledge or new to the programs.
  * Programs should use multiple threads whenever possible or feasible.
  * C# programs should be fully compatible with both:
    * Mono 2.10.x or later
    * The latest release of the .NET Framework
  * SCons is used to build programs (<font color='#ff0000'>experimental</font>)

## Dependencies ##
  * C++11 or C# 4.0 compatible compiler (with corresponding framework version)
  * libxml2 (and its dependencies)
  * zlib >= 1.28
  * bzip
  * Eigen >= 3.2.x
  * Boost >= 1.54.x
  * SCons >= 2.3.0

## Under consideration ##
  * Replace zlib with liblzma.
  * Write binary instead of plain text.
  * Use SQLite or BerkeleyDB for trajectories
  * One GUI for all programs
  * Compatibility with C++14 (<font color='#ff0000'>as soon as compilers are available</font>)
  * Yeppp!
  * HDF5
  * Boost.SIMD
  * Boost.Compute

## General stuff on to-do list ##
  * Implement a thread pool using features introduced by C++11
  * Test Intel® SPMD Program Compiler 1.6.0
  * Testing HDF5