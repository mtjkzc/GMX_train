#!/bin/bash

###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################
#	███╗   ███╗██████╗     ███████╗ ██████╗██████╗ ██╗██████╗ ████████╗
#	████╗ ████║██╔══██╗    ██╔════╝██╔════╝██╔══██╗██║██╔══██╗╚══██╔══╝
#	██╔████╔██║██║  ██║    ███████╗██║     ██████╔╝██║██████╔╝   ██║   
#	██║╚██╔╝██║██║  ██║    ╚════██║██║     ██╔══██╗██║██╔═══╝    ██║   
#	██║ ╚═╝ ██║██████╔╝    ███████║╚██████╗██║  ██║██║██║        ██║   
#	╚═╝     ╚═╝╚═════╝     ╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝        ╚═╝  
###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################

###  MD  pipeline  main  script  ###

# v0.2.1 - DNA Project
# 2024.06.13
# University of Zagreb, Facoulty of Science, Department of Chemistry
# Matej Kožić | mkozic@chem.pmf.hr

# master.sh runs this script to run MD simulations
# anls_script.sh follows after this script to analyze iz
# graph_report.py creates a graphical report after it

# Workflow:
#	Step 0: Adjusting values, parameters, & variables
#	Step 1: System solvatation and topology preparation
#	Step 2: Neutralization, minimization, NVT equilibration, & NPT equilibration
#	Step 3: MD production
#	END

# End of comment

###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###############################################################################
#	███████╗████████╗███████╗██████╗      ██████╗ 
#	██╔════╝╚══██╔══╝██╔════╝██╔══██╗    ██╔═████╗
#	███████╗   ██║   █████╗  ██████╔╝    ██║██╔██║
#	╚════██║   ██║   ██╔══╝  ██╔═══╝     ████╔╝██║
#	███████║   ██║   ███████╗██║         ╚██████╔╝
#	╚══════╝   ╚═╝   ╚══════╝╚═╝          ╚═════╝ 
###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###############################################################################

###--------------------------------------------------###
###  Step  0:  Assigning  variables  &  MDP  values  ###
###--------------------------------------------------###

### Setup variables ###

# Run name
run_name=$1

# Forcefield selection
ff=1  	# Will select forcefield from the dirtectory when 1

# Setup  variables
boxsize=1
ignore_hydrogen="yes"
keep_crys_wat="no"

#  System  type
number_of_dna_chains=2

# Settings  for  mdrun
gpuid=0


###  MDP  values  ###

# Minimization #
#Number of minimization steps
min_nsteps=1000

# NVT equilibration #
#NVT timestep (dt), number of steps, writing frequency
nvt_dt=0.001 #1fs
nvt_nsteps=600000
nvt_nstxtcout=5000
nvt_nstcalcenergy=100
nvt_nstenergy=5000 
nvt_nstlog=5000

#NVT temperature coupling
nvt_tc_group1="non-Water"
nvt_tc_group2="Water"
nvt_ref_temp=300
nvt_annealing_time=500
nvt_annealing_temp_start=10
nvt_annealing_temp_goal=300

# NPT equilibration #
#NPT timestep (dt), number of steps, writing frequency
npt_dt=0.001 #1fs
npt_nsteps=500000
npt_nstxtcout=5000
npt_nstcalcenergy=100
npt_nstenergy=5000
npt_nstlog=5000

#NPT temperature coupling
npt_tc_group1="non-Water"
npt_tc_group2="Water"
npt_ref_temp=300

# NPT MD production #
#NPT MD timestep, number of steps, write freq, temp
npt_md_dt=0.002 #2fs
npt_md_nsteps=2500000	#250000 #5 ns; 250 000 000 #500ns
npt_md_nstxtcout=5000
npt_md_nstcalcenerg=100
npt_md_nstenergy=1000
npt_md_nstlog=1000
npt_md_ref_temp=300


### Settings for trjconv ###

trjconv_skip=1


### Master script variables ###

name="-"
logfile=${name}.log

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###  Finding  a  PDB  file  in  current  directory  ###

directory="."
extension=".pdb" 

if [ -d "$directory" ]; then
    for file in "$directory"/*"$extension"; do
        if [ -e "$file" ]; then
            name="$file"
            break
        fi
    done
    if [ -n "$name" ]; then
        echo "Found file with extension '$extension': $name"
    else
        echo "No files with extension '$extension' found in $directory"
    fi
else
    echo "Directory '$directory' does not exist."
fi

echo "= = = PDB file: "$name" = = ="

name=${name/.pdb/}


###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
#	███████╗████████╗███████╗██████╗      ██╗
#	██╔════╝╚══██╔══╝██╔════╝██╔══██╗    ███║
#	███████╗   ██║   █████╗  ██████╔╝    ╚██║
#	╚════██║   ██║   ██╔══╝  ██╔═══╝      ██║
#	███████║   ██║   ███████╗██║          ██║
#	╚══════╝   ╚═╝   ╚══════╝╚═╝          ╚═╝
###############################################################################                                                                              
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###-------------------------------------------###
###  Step  1:  System  topology  preparation  ###
###-------------------------------------------###

###  Creating  an  index  file  ###

echo "|------------------------------------------------|"
echo "|= = = = = = = = = = INDEXING = = = = = = = = = =|"
echo "|------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

gmx pdb2gmx -f ${name}.pdb -o ${name}_HYDROGENIZED.pdb -p ${name}_temp.top -n ${name}_CHAINS.ndx -water tip3p << EOF
$ff
EOF

gmx make_ndx -f ${name}_HYDROGENIZED.pdb -o ${name}_CHAINS.ndx << EOF
chain A
chain B
q
EOF


###  Topology  &  solvataion  ###

echo "|---------------------------------------------------|"
echo "|= = = = = = = = = = SOLVATATION = = = = = = = = = =|"
echo "|---------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

gmx pdb2gmx -f ${name}_HYDROGENIZED.pdb -o ${name}.gro -p ${name}.top -n ${name}.ndx -water tip3p << EOF
$ff
EOF

gmx editconf -f ${name}.gro -bt cubic -d $boxsize -o ${name}_box.gro

gmx solvate -cp ${name}_box.gro -cs spc216.gro -p ${name}.top -o ${name}_sol.gro

mkdir prep
cp ${name}.gro ${name}.top ${name}_box.gro ${name}_sol.gro prep

echo "= = = Finished solvatation for: "${name}_sol.gro" = = ="


###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
#	███████╗████████╗███████╗██████╗     ██████╗
#	██╔════╝╚══██╔══╝██╔════╝██╔══██╗    ╚════██╗
#	███████╗   ██║   █████╗  ██████╔╝     █████╔╝
#	╚════██║   ██║   ██╔══╝  ██╔═══╝     ██╔═══╝
#	███████║   ██║   ███████╗██║         ███████╗
#	╚══════╝   ╚═╝   ╚══════╝╚═╝         ╚══════╝
###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###--------------------------------------------------------------###
###  Step  2:  Neutralization,  Minimization,  &  equilibration  ###
###--------------------------------------------------------------###

###  Modifying  MDP  scripts  ###

sed -i "s/min_nsteps/${min_nsteps}/g" mdp_py_sh/min.mdp

sed -i "s/nvt_dt/${nvt_dt}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_nsteps/${nvt_nsteps}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_nstxtcout/${nvt_nstxtcout}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_nstcalcenergy/${nvt_nstcalcenergy}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_nstenergy/${nvt_nstenergy}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_nstlog/${nvt_nstlog}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_tc_group1/${nvt_tc_group1}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_tc_group2/${nvt_tc_group2}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_ref_temp/${nvt_ref_temp}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_annealing_time/${nvt_annealing_time}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_annealing_temp_start/${nvt_annealing_temp_start}/g" mdp_py_sh/nvt.mdp
sed -i "s/nvt_annealing_temp_goal/${nvt_annealing_temp_goal}/g" mdp_py_sh/nvt.mdp

sed -i "s/npt_dt/${npt_dt}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_nsteps/${npt_nsteps}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_nstxtcout/${npt_nstxtcout}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_nstcalcenergy/${npt_nstcalcenergy}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_nstenergy/${npt_nstenergy}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_nstlog/${npt_nstlog}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_tc_group1/${npt_tc_group1}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_tc_group2/${npt_tc_group2}/g" mdp_py_sh/npt.mdp
sed -i "s/npt_ref_temp/${npt_ref_temp}/g" mdp_py_sh/npt.mdp

sed -i "s/npt_md_dt/${npt_md_dt}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_nsteps/${npt_md_nsteps}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_nstxtcout/${npt_md_nstxtcout}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_nstcalcenerg/${npt_md_nstcalcenerg}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_nstenergy/${npt_md_nstenergy}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_nstlog/${npt_md_nstlog}/g" mdp_py_sh/npt_md.mdp
sed -i "s/npt_md_ref_temp/${npt_md_ref_temp}/g" mdp_py_sh/npt_md.mdp

sed -i "s/annealing_time_end/${annealing_time_end}/g" mdp_py_sh/npt_md.mdp
sed -i "s/annealing_temp_start/${annealing_temp_start}/g" mdp_py_sh/npt_md.mdp
sed -i "s/annealing_temp_goal/${annealing_temp_goal}/g" mdp_py_sh/npt_md.mdp


###  Modifying  analysis  script  ###

sed -i "s/_trjconv_skip_/${trjconv_skip}/g" ./anls_script.sh

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###  Neutralization  ###

echo "|------------------------------------------------------|"
echo "|= = = = = = = = = = NEUTRALIZATION = = = = = = = = = =|"
echo "|------------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

mkdir analysis neut min nvt npt md

mv ./mdp_py_sh/min.mdp .

gmx grompp -f min.mdp -c ${name}_sol.gro -r ${name}_sol.gro -p ${name}.top -o ${name}_ion.tpr -maxwarn 1

echo 3 | gmx genion -p ${name}.top -s ${name}_ion.tpr -o ${name}_ion.gro -pname NA -nname CL -neutral

mv ${name}_sol.gro ${name}_ion.tpr neut

echo "= = = Finished neutralization = = ="

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###  Minimization  ###

echo "|----------------------------------------------------|"
echo "|= = = = = = = = = = MINIMIZATION = = = = = = = = = =|"
echo "|----------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

gmx grompp -f min.mdp -c ${name}_ion.gro -r ${name}_ion.gro -p ${name}.top -o ${name}_min.tpr 

gmx mdrun -v -s ${name}_min.tpr -o ${name}_min.trr -c ${name}_min.gro -e ${name}_min.edr -gpu_id ${gpuid} > min_output.log 2>&1

tail min_output.log


echo 10 | gmx energy -f ${name}_min.edr -o analysis/min_potential.xvg

mv ${name}_ion.gro ${name}_min.tpr min.mdp ${name}_min.trr ${name}_min.edr min

echo "= = = Finished minimization = = ="

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###  NVT  equilibration  ###

echo "|---------------------------------------------------------|"
echo "|= = = = = = = = = = NVT equilibration = = = = = = = = = =|"
echo "|---------------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

mv ./mdp_py_sh/nvt.mdp .

gmx grompp -f nvt.mdp -c ${name}_min.gro -r ${name}_min.gro -p ${name}.top -o ${name}_nvt.tpr

gmx mdrun -nb gpu -bonded gpu -pme gpu -update gpu -v -s ${name}_nvt.tpr -o ${name}_nvt.trr -c ${name}_nvt.gro -e ${name}_nvt.edr -gpu_id ${gpuid} > nvt_output.log 2>&1

tail nvt_output.log

echo 14 | gmx energy -f ${name}_nvt.edr -o analysis/nvt_temperature.xvg
echo 15 | gmx energy -f ${name}_nvt.edr -o analysis/nvt_pressure.xvg
echo 10 | gmx energy -f ${name}_nvt.edr -o analysis/nvt_potential.xvg

mv nvt.mdp ${name}_min.gro ${name}_nvt.tpr ${name}_nvt.edr ${name}_nvt.trr nvt

echo "= = = Finished NVT equilibration = = ="

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

###  NPT  equilibration  ###

echo "|---------------------------------------------------------|"
echo "|= = = = = = = = = = NPT equilibration = = = = = = = = = =|"
echo "|---------------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

mv ./mdp_py_sh/npt.mdp .

gmx grompp -f npt.mdp -c ${name}_nvt.gro -r ${name}_nvt.gro -p ${name}.top -o ${name}_npt.tpr

gmx mdrun -nb gpu -bonded gpu -pme gpu -update gpu -v -s ${name}_npt.tpr -o ${name}_npt.trr -x ${name}_rawtraj.xtc -c ${name}_npt.gro -e ${name}_npt.edr -gpu_id ${gpuid} > npt_output.log 2>&1

tail npt_output.log

echo 14 | gmx energy -f ${name}_npt.edr -o analysis/npt_temperature.xvg
echo 15 | gmx energy -f ${name}_npt.edr -o analysis/npt_pressure.xvg
echo 10 | gmx energy -f ${name}_npt.edr -o analysis/npt_potential.xvg

mv npt.mdp ${name}_nvt.gro ${name}_npt.tpr ${name}_npt.trr npt

echo "= = = Finished NPT equilibration = = ="

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
#	███████╗████████╗███████╗██████╗     ██████╗ 
#	██╔════╝╚══██╔══╝██╔════╝██╔══██╗    ╚════██╗
#	███████╗   ██║   █████╗  ██████╔╝     █████╔╝
#	╚════██║   ██║   ██╔══╝  ██╔═══╝      ╚═══██╗
#	███████║   ██║   ███████╗██║         ██████╔╝
#	╚══════╝   ╚═╝   ╚══════╝╚═╝         ╚═════╝ 
###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################


###----------------------------###
###  Step  4:  MD  production  ###
###----------------------------###

echo "|-----------------------------------------------------|"
echo "|= = = = = = = = = = MD PRODUCTION = = = = = = = = = =|"
echo "|-----------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

mkdir dump
mv *.xtc* *.cpt* dump

mv ./mdp_py_sh/npt_md.mdp .

gmx grompp -f npt_md.mdp -c ${name}_npt.gro -r ${name}_npt.gro -p ${name}.top -o ${name}_md.tpr

gmx mdrun -nb gpu -bonded gpu -pme gpu -update gpu -v -s ${name}_md.tpr -o ${name}_md.trr -e ${name}_md.edr -x ${name}_rawtraj.xtc -gpu_id ${gpuid} > md_output.log 2>&1

tail md_output.log

echo "= = = Finished MD production = = ="

echo "|---------------------------------------|"
echo "| = = = = = = = =       = = = = = = = = |"
echo "| > > > > > > MD SCRIPT END < < < < < < |"
echo "| = = = = = = = =       = = = = = = = = |"
echo "|---------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

###############################################################################
# END END END END END END END END END END END END END END END END END END END #
###############################################################################
#	███████╗███╗   ██╗██████╗ 
#	██╔════╝████╗  ██║██╔══██╗
#	█████╗  ██╔██╗ ██║██║  ██║
#	██╔══╝  ██║╚██╗██║██║  ██║
#	███████╗██║ ╚████║██████╔╝
#	╚══════╝╚═╝  ╚═══╝╚═════╝
###############################################################################
# END END END END END END END END END END END END END END END END END END END #
### E N D #####################################################################

