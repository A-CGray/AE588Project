#!/bin/bash

mkdir -p scalingTest

for alg in  SNOPT IPOPT NLPQLP #ParOptSL1 ParOptFilter ParOptMMA
do
    
    mkdir -p scalingTest/"$alg"
    if ["$alg" = "SNOPT"];then
        hess=50
    else
        hess=1000
    fi
    rm -rf scalingTest/"$alg"/4
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 1 --nSpar 1 --nRib 1 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/0 --optimiser "$alg"
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 9 --nSpar 9 --nRib 9 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/1 --optimiser "$alg"
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 18 --nSpar 18 --nRib 18 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/2 --optimiser "$alg"
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 36 --nSpar 36 --nRib 18 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/3 --optimiser "$alg"
    mpirun -np 4 python struct_opt.py --task opt --nSkin 72 --nSpar 36 --nRib 18 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/4 --optimiser "$alg" --saveOptIters major
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 180 --nSpar 36 --nRib 18 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/5 --optimiser "$alg"
    # mpirun -np 4 python struct_opt.py --task opt --nSkin 360 --nSpar 36 --nRib 18 --adjCon 2.5e-3 --hessianUpdate "$hess"  --optIter 600 --output scalingTest/"$alg"/6 --optimiser "$alg" --saveOptIters major
done