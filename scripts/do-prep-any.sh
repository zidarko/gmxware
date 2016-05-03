#!/bin/bash

# This script will construct and prepare a ligand-protein system for simulation
# NB: 
# * input molecule must be in ./input-mols/
# * the following files are needed: molecule.{prm,itp}, molecule_ini.pdb
# (the mandatory files can be generated using cgenff_charmm2gmx.py script)
# * supplementary files will be copied over for consistency: molecule.{mol2,str}

# Eh ...
echo Enter name of the molecule.
read input

# Less changes to be made.
i="$input"

# 1. Copy ligand parameter data
	mkdir "$i"
	cp -av ./input-mols/"$i".mol2 ./"$i"/
	cp -av ./input-mols/"$i".str ./"$i"/
	cp -av ./input-mols/"$i".prm ./"$i"/
	cp -av ./input-mols/"$i".itp ./"$i"/
	cp -av ./input-mols/"$i"_ini.pdb ./"$i"/
# 2. Generate appropriate GRO file for the receptor
	cp -av ./template/charmm* ./"$i"/
	cp -av ./template/er-a-a.itp ./"$i"/
	cp -av ./template/er-a-a-posres.itp ./"$i"/
	head -n4005 ./template/er-a-a.gro > ./"$i"/temp-rec.gro
# 3. Generate appropriate GRO for ligand and merge it with receptor file
	gmx genconf -f ./"$i"/"$i"_ini.pdb -o ./"$i"/"$i"_ini.gro
	echo 1 | gmx genrestr -f ./"$i"/"$i"_ini.gro -o ./"$i"/"$i"-posres
	sed '$d' < ./"$i"/"$i"_ini.gro > ./"$i"/temp-lig.gro
	sed -i -e 1,2d ./"$i"/temp-lig.gro
# 3b. Adjust number of atoms in the future receptor-ligand GRO file
	lign=`wc ./"$i"/temp-lig.gro | awk '{print $1}'`
	recn=`head -n2 ./"$i"/temp-rec.gro | tail -n1`
	totn=`echo "$recn + $lign"| bc`
	sed -i "s/$recn/$totn/" ./"$i"/temp-rec.gro
	cat ./"$i"/temp-rec.gro ./"$i"/temp-lig.gro > ./"$i"/er-a-a-"$i".gro
	tail -n1 ./template/er-a-a.gro >> ./"$i"/er-a-a-"$i".gro
	gmx genconf -f ./"$i"/er-a-a-"$i".gro -o ./"$i"/er-a-a-"$i".gro -renumber
	var=`echo "$i" | tr [a-z] [A-Z]`
	sed -i "s/template/$var/" ./"$i"/er-a-a-"$i".gro
# 3c. Change ligand name to LIG (makes life easier later on)
	cp -av ./"$i"/"$i".itp ./"$i"/esdi.itp.bak
	sed -i "s/$var/LIG/" ./"$i"/"$i".itp
# 4. Generate appropriate TOP file using a topo-main-template.top
	sed "s/template/$i/" < ./template/topo-main-template.top > ./"$i"/er-a-a-"$i".top
	sed -i "s/temp2/$var/" ./"$i"/er-a-a-"$i".top
# 5. Convert initial complex to PDB for visual inspection
	gmx genconf -f ./"$i"/er-a-a-"$i".gro -o ./"$i"/er-a-a-"$i".pdb
# 6. Adjust the box size, solvate and add ions
	gmx editconf -f ./"$i"/er-a-a-"$i".gro -o ./"$i"/er-a-a-"$i".gro -bt cubic -d 0.5 -c
	cp -av ./"$i"/er-a-a-"$i".top ./"$i"/wer-a-a-"$i".top
	gmx solvate 	-cp ./"$i"/er-a-a-"$i".gro -cs spc216.gro \
			-p ./"$i"/wer-a-a-"$i".top \
			-o ./"$i"/wer-a-a-"$i".gro
	gmx grompp	-f ./input/em-fix.mdp \
			-c ./"$i"/wer-a-a-"$i"\
			-p ./"$i"/wer-a-a-"$i"\
			-o ./"$i"/temp -po ./"$i"/temp
	cp -av ./"$i"/wer-a-a-"$i".top ./"$i"/wer-a-a-"$i"-ion.top
	echo 15 |	gmx genion 	-s ./"$i"/temp \
					-p ./"$i"/wer-a-a-"$i"-ion.top \
					-o ./"$i"/wer-a-a-"$i"-ion.gro \
					-neutral -conc 0.154
# 7. Generate index file with groups neccessary for MM/PBSA
	gmx make_ndx 	-f ./"$i"/wer-a-a-"$i"-ion.gro \
			-o ./"$i"/wer-a-a-"$i"-ion <<EOF
1|13
13
name 25 blah
del 20
q
EOF
	sed -i "s/blah/$var/" ./"$i"/wer-a-a-"$i"-ion.ndx
# 8. Compile EM run file
	gmx grompp	-f ./input/em-fix.mdp \
			-c ./"$i"/wer-a-a-"$i"-ion \
			-p ./"$i"/wer-a-a-"$i"-ion \
			-n ./"$i"/wer-a-a-"$i"-ion \
			-o ./"$i"/wer-a-a-"$i"-ion-em-fix -po ./"$i"/temp
# 9. Clean-up
	rm ./"$i"/temp* ./"$i"/\#*
	rm \#*
# 10. Check if intial TPR files is present, print error if not.
	init=wer-a-a-"$i"-ion-em-fix.tpr
	echo "${init/.tpr}": `if [ ! -e ./"$i"/"$init" ] ; then echo -e Not Found ; else echo -e Found ; fi`



