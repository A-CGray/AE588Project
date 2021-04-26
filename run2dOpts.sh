#!/bin/bash

mkdir -p 2dOpts

for alg in SNOPT IPOPT NLPQLP ParOptSL1 ParOptFilter ParOptMMA
do
    mpirun -np 4 python struct_opt.py --task 2dOpt --output 2dOpts/"$alg" --optimiser "$alg"
done