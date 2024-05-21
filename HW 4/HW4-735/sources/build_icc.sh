#!/bin/bash
module load intel
mpiicpc -o qsort_hypercube.exe qsort_hypercube.cpp
mpiicpc -o qsort_hypercube_descending.exe qsort_hypercube_descending.cpp
