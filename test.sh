#!/bin/bash

mpicc -o lab06.exe lab06.c -lm
mpiexec -np 1 ./lab06.exe
