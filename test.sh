#!/bin/bash

mpicc -o finals.exe finals.c -lm
mpiexec -np 1 ./finals.exe
mpiexec -np 1 ./finals.exe
mpiexec -np 1 ./finals.exe
mpiexec -np 2 ./finals.exe
mpiexec -np 2 ./finals.exe
mpiexec -np 2 ./finals.exe
mpiexec -np 3 ./finals.exe
mpiexec -np 3 ./finals.exe
mpiexec -np 3 ./finals.exe
mpiexec -np 4 ./finals.exe
mpiexec -np 4 ./finals.exe
mpiexec -np 4 ./finals.exe