#!/bin/bash

mkdir -p MultiModalTest
rm -rf MultiModalTest/*

# Run multimodality test with and without adjacency constraints
mkdir -p MultiModalTest/AdjCon
mkdir -p MultiModalTest/NoAdjCon

# First run at uniform initial design
mpirun -np 4 python struct_opt.py --output MultiModalTest/AdjCon/0 --task opt --nSkin 72 --nSpar 36 --nRib 18 --adjCon 5e-3 --optimiser SNOPT --hessianUpdate 50 --saveOptIters major --optIter 600

mpirun -np 4 python struct_opt.py --output MultiModalTest/NoAdjCon/0 --task opt --nSkin 72 --nSpar 36 --nRib 18 --optimiser SNOPT --hessianUpdate 50 --saveOptIters major --optIter 600

# Now run from 4 random starting points
for i in {1..4}
	do

	if (( $i == 1 )); then
		save="major"
	else
		save="none"
	fi
	    mpirun -np 4 python struct_opt.py --output MultiModalTest/AdjCon/"$i" --task opt --nSkin 72 --nSpar 36 --nRib 18 --adjCon 5e-3 --optimiser SNOPT --randomX0 --hessianUpdate 50 --saveOptIters "$save" --optIter 600

	    mpirun -np 4 python struct_opt.py --output MultiModalTest/NoAdjCon/"$i" --task opt --nSkin 72 --nSpar 36 --nRib 18 --optimiser SNOPT --randomX0 --hessianUpdate 50 --saveOptIters "$save" --optIter 600
done