# simple-adaptive-sparse-grids

A single header C++ implementation of multi-dimensional interpolating sparse grids, featuring adaptive refinement and MATLAB support.

![Alt Text](https://github.com/ArthurFirmino/simple-adaptive-sparse-grids/raw/master/example/figure.png)

## Background

This implementation is based in part on the present literature surrounding sparse grids, with its main motivation being simple and 
concise code that is easy to understand. 

**simple-adaptive-sparse-grids** doesn't feature any specialized algorithms on sparse grids, 
having only hierarchization, evaluation, refinement, and unrefinment operations.

Grid nodes are stored in a hash table, similar to the [SG++](https://github.com/SGpp/SGpp) library which does implement
other algorithms on sparse grids.
