#!/bin/bash

# script to run equilibration of protein-ligand system
for i in molecule
	do
	input=wer-a-a-"$i"-ion
	for curr in nvt-fix npt-fix npt-prod
		do
		if [ "$curr" == nvt-fix ]
			then
			prev=em-fix
		fi
		if [ "$curr" == npt-fix ]
			then
			prev=nvt-fix
		fi
		if [ "$curr" == npt-prod ]
			then
			prev=npt-fix
		fi
		cd ./"$i"
		gmx grompp	-f ../input/"$curr".mdp \
				-c "$input" \
				-p "$input" \
				-n "$input" \
				-o "$input"-"$curr" \
				-t "$input"-"$prev".trr \
				-po temp
		rm temp.mdp
		if [ "$curr" == npt-prod ]
			then
			gmx mdrun	-v -deffnm "$input"-"$curr" \
					-pin on -ntomp 5 -gpu_id 0011 \
					-nsteps 25000000
			else
			gmx mdrun	-v -deffnm "$input"-"$curr" \
					-pin on -ntomp 5 -gpu_id 0011
		fi
		rm \#*
		cd ../
		done
	done
