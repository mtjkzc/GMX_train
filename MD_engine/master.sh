#!/bin/bash

###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################
#	███╗   ███╗ █████╗ ███████╗████████╗███████╗██████╗ 
#	████╗ ████║██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
#	██╔████╔██║███████║███████╗   ██║   █████╗  ██████╔╝
#	██║╚██╔╝██║██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
#	██║ ╚═╝ ██║██║  ██║███████║   ██║   ███████╗██║  ██║
#	╚═╝     ╚═╝╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝
###############################################################################                                                  
# START START START START START START START START START START START START STA #
###############################################################################

# 2024.06.12 
# University of Zagreb, Facoulty of Science, Department of Chemistry
# Matej Kožić | mkozic@chem.pmf.hr 

# A script for running MD simulations starting from AlphaFold 3 structures.
# It runs 3 replicas, starting from 3 models (m0, m1, m2) from AlphaFold.
# AlphaFold structures need to be named "name"_af_m0.cif, ...m2.cif .... 

# Requirements to be in a dir same as this script:
#	- AlphaFold structures
#	- master.sh (This script)
#	- md_script.sh
#	- mdp_py_sh directory

# End of comment

###############################################################################
#00000000000000000000000000000000000000000000000000000000000000000000000000000#
###############################################################################

###------------------------------###
###  Step  0:  Conversion  Hell  ###
###------------------------------###

### Assigning name ###
core_name=$1
name="${core_name}_af"

### Naming directories ###

dir0="replica_0"
dir1="replica_1"
dir2="replica_2"


##--------------------------##
###  Processing  model  0  ###
##--------------------------##

# Extracting AlphaFold models from .zip
mv mdp_py_sh/af_extractor.sh .
./af_extractor.sh "$core_name"


# Removing 5' terminal phosphate cap
mv mdp_py_sh/dephosphate.sh .
./dephosphate.sh "./${core_name}_af_m0_PHOSPHATED.cif"
./dephosphate.sh "./${core_name}_af_m1_PHOSPHATED.cif"
./dephosphate.sh "./${core_name}_af_m2_PHOSPHATED.cif"

mv ${core_name}_af_m0_PHOSPHATED_CLEANED.cif ${name}_m0.cif
mv ${core_name}_af_m1_PHOSPHATED_CLEANED.cif ${name}_m1.cif
mv ${core_name}_af_m2_PHOSPHATED_CLEANED.cif ${name}_m2.cif


# Converting AlphaFold .cif to AlphaFold noncharmm .pdb
cif_filename=${name}_m0.cif
pdb_savename=${name}_m0_noncharmm.pdb

cp mdp_py_sh/cif2pdb.py ./cif2pdb_temp.py
sed -i "s/cif_filename/${cif_filename}/g" cif2pdb_temp.py
sed -i "s/output_pdb_filename/${pdb_savename}/g" cif2pdb_temp.py

python3 ./cif2pdb_temp.py
rm ./cif2pdb_temp.py

# Converting AlphaFold PDB to CHARMM nomenclature
input_file=${name}_m0_noncharmm.pdb
output_file=${name}_m0.pdb

# Create a temporary file for intermediate processing
temp_file=$(mktemp)

# Perform the renaming operations
sed -e 's/\bOP1\b/O1P/g' \
    -e 's/\bOP2\b/O2P/g' \
    -e '/ OP3 /d' \
    "$input_file" > "$temp_file"
    
# Move the temp file to the output file
mv "$temp_file" "$output_file"

echo "========== Renaming completed. Written to $output_file =========="


##--------------------------##
###  Processing  model  1  ###
##--------------------------##

# Converting AlphaFold .cif to AlphaFold noncharmm .pdb
cif_filename=${name}_m1.cif
pdb_savename=${name}_m1_noncharmm.pdb

cp mdp_py_sh/cif2pdb.py ./cif2pdb_temp.py
sed -i "s/cif_filename/${cif_filename}/g" cif2pdb_temp.py
sed -i "s/output_pdb_filename/${pdb_savename}/g" cif2pdb_temp.py

python3 ./cif2pdb_temp.py
rm ./cif2pdb_temp.py

# Converting AlphaFold PDB to CHARMM nomenclature
input_file=${name}_m1_noncharmm.pdb
output_file=${name}_m1.pdb

# Create a temporary file for intermediate processing
temp_file=$(mktemp)

# Perform the renaming operations
sed -e 's/\bOP1\b/O1P/g' \
    -e 's/\bOP2\b/O2P/g' \
    -e '/ OP3 /d' \
    "$input_file" > "$temp_file"
    
# Move the temp file to the output file
mv "$temp_file" "$output_file"

echo "========== Renaming completed. Written to $output_file =========="


##--------------------------##
###  Processing  model  2  ###
##--------------------------##

# Converting AlphaFold .cif to AlphaFold noncharmm .pdb
cif_filename=${name}_m2.cif
pdb_savename=${name}_m2_noncharmm.pdb

cp mdp_py_sh/cif2pdb.py ./cif2pdb_temp.py
sed -i "s/cif_filename/${cif_filename}/g" cif2pdb_temp.py
sed -i "s/output_pdb_filename/${pdb_savename}/g" cif2pdb_temp.py

python3 ./cif2pdb_temp.py
rm ./cif2pdb_temp.py

# Converting AlphaFold PDB to CHARMM nomenclature
input_file=${name}_m2_noncharmm.pdb
output_file=${name}_m2.pdb

# Create a temporary file for intermediate processing
temp_file=$(mktemp)

# Perform the renaming operations
sed -e 's/\bOP1\b/O1P/g' \
    -e 's/\bOP2\b/O2P/g' \
    -e '/ OP3 /d' \
    "$input_file" > "$temp_file"
    
# Move the temp file to the output file
mv "$temp_file" "$output_file"

echo "========== Renaming completed. Written to $output_file =========="


###############################################################################
#11111111111111111111111111111111111111111111111111111111111111111111111111111#
###############################################################################

###---------------------------------------###
###  Step  1:  Preparation  purgatorium   ###
###---------------------------------------###


### Turning on Gromacs ###

source /opt/gromacs/bin/GMXRC.bash


### Creating Directories ###

mkdir $dir0
mkdir $dir1
mkdir $dir2

mkdir ${name}_MD_output
mkdir ${name}_MD_graphs


### Moving things into MD directories ###

# Forcefields
cp -r *.ff $dir0
cp -r *.ff $dir1
mv *.ff $dir2

# Structures
mv ${name}_m0.pdb $dir0
mv ${name}_m1.pdb $dir1
mv ${name}_m2.pdb $dir2

# MD scripts
cp md_script.sh $dir0
cp md_script.sh $dir1
mv md_script.sh $dir2

# Analysis script
cp mdp_py_sh/anls_script.sh $dir0
cp mdp_py_sh/anls_script.sh $dir1
mv mdp_py_sh/anls_script.sh $dir2

# Graph report script
mv mdp_py_sh/graph_report.py .

# MDP directory
cp -r mdp_py_sh $dir0
cp -r mdp_py_sh $dir1
mv mdp_py_sh $dir2


###############################################################################
#22222222222222222222222222222222222222222222222222222222222222222222222222222#
###############################################################################

###---------------------------------------------------###
###  Step  2:  MD  Initialization  circle  of  hell   ###
###---------------------------------------------------###


### MD of replica 0 from model 0 ###

echo "|------------------------------------------------------|"
echo "|= = = = = = = = = =  MD  MODEL  0  = = = = = = = = = =|" 
echo "|------------------------------------------------------|"

(cd $dir0 && ./md_script.sh "$name | REPLICA 0")
(cd $dir0 && ./anls_script.sh "$name | REPLICA 0")
cp -r ${dir0}/output_md ./output_${dir0}
cp ${dir0}/analysis/* ./output_${dir0}


### MD of replica 1 from model 1 ###

echo "|------------------------------------------------------|"
echo "|= = = = = = = = = =  MD  MODEL  1  = = = = = = = = = =|" 
echo "|------------------------------------------------------|"

(cd $dir1 && ./md_script.sh "$name | REPLICA 1")
(cd $dir1 && ./anls_script.sh "$name | REPLICA 1")
cp -r ${dir1}/output_md ./output_${dir1}
cp ${dir1}/analysis/* ./output_${dir1}


### MD of replica 1 from model 1 ###

echo "|------------------------------------------------------|"
echo "|= = = = = = = = = =  MD  MODEL  2  = = = = = = = = = =|" 
echo "|------------------------------------------------------|"

(cd $dir2 && ./md_script.sh "$name | REPLICA 2")
(cd $dir2 && ./anls_script.sh "$name | REPLICA 2")
cp -r ${dir2}/output_md ./output_${dir2}
cp ${dir2}/analysis/* ./output_${dir2}


###############################################################################
#33333333333333333333333333333333333333333333333333333333333333333333333333333#
###############################################################################

###------------------------------###
###  Step  3:  Analysis   heaven ###
###------------------------------###


### Rounding up all MD trajectories ###

# Trajectories
cp output_${dir0}/*.xtc ./${name}_MD_output
cp output_${dir1}/*.xtc ./${name}_MD_output
cp output_${dir2}/*.xtc ./${name}_MD_output

# Topologies
cp output_${dir0}/*.gro ./${name}_MD_output
cp output_${dir1}/*.gro ./${name}_MD_output
cp output_${dir2}/*.gro ./${name}_MD_output

echo "========== Moved all trajectories into MD_output =========="


### Rounding up all graphs ###

a1="min_potential.xvg"
b1="nvt_potential.xvg"
b2="nvt_pressure.xvg"
b3="nvt_temperature.xvg"
c1="npt_potential.xvg"
c2="npt_pressure.xvg"
c3="npt_temperature.xvg"
d1="md_potential.xvg"
d2="md_pressure.xvg"
d3="md_temperature.xvg"
e1="rmsd.xvg"
e2="rmsf.xvg"
e3="gyr.xvg"


# Replica 0 with model 0 
cp output_${dir0}/${a1} ./${name}_MD_graphs/${name}_m0_${a1}
cp output_${dir0}/${b1} ./${name}_MD_graphs/${name}_m0_${b1}
cp output_${dir0}/${b2} ./${name}_MD_graphs/${name}_m0_${b2}
cp output_${dir0}/${b3} ./${name}_MD_graphs/${name}_m0_${b3}
cp output_${dir0}/${c1} ./${name}_MD_graphs/${name}_m0_${c1}
cp output_${dir0}/${c2} ./${name}_MD_graphs/${name}_m0_${c2}
cp output_${dir0}/${c3} ./${name}_MD_graphs/${name}_m0_${c3}
cp output_${dir0}/${d1} ./${name}_MD_graphs/${name}_m0_${d1}
cp output_${dir0}/${d2} ./${name}_MD_graphs/${name}_m0_${d2}
cp output_${dir0}/${d3} ./${name}_MD_graphs/${name}_m0_${d3}
cp output_${dir0}/${name}_m0_${e1} ./${name}_MD_graphs
cp output_${dir0}/${name}_m0_${e2} ./${name}_MD_graphs
cp output_${dir0}/${name}_m0_${e3} ./${name}_MD_graphs


# Replica 1 with model 1
cp output_${dir1}/${a1} ./${name}_MD_graphs/${name}_m1_${a1}
cp output_${dir1}/${b1} ./${name}_MD_graphs/${name}_m1_${b1}
cp output_${dir1}/${b2} ./${name}_MD_graphs/${name}_m1_${b2}
cp output_${dir1}/${b3} ./${name}_MD_graphs/${name}_m1_${b3}
cp output_${dir1}/${c1} ./${name}_MD_graphs/${name}_m1_${c1}
cp output_${dir1}/${c2} ./${name}_MD_graphs/${name}_m1_${c2}
cp output_${dir1}/${c3} ./${name}_MD_graphs/${name}_m1_${c3}
cp output_${dir1}/${d1} ./${name}_MD_graphs/${name}_m1_${d1}
cp output_${dir1}/${d2} ./${name}_MD_graphs/${name}_m1_${d2}
cp output_${dir1}/${d3} ./${name}_MD_graphs/${name}_m1_${d3}
cp output_${dir1}/${name}_m1_${e1} ./${name}_MD_graphs
cp output_${dir1}/${name}_m1_${e2} ./${name}_MD_graphs
cp output_${dir1}/${name}_m1_${e3} ./${name}_MD_graphs


# Replica 2 with model 2
cp output_${dir2}/${a1} ./${name}_MD_graphs/${name}_m2_${a1}
cp output_${dir2}/${b1} ./${name}_MD_graphs/${name}_m2_${b1}
cp output_${dir2}/${b2} ./${name}_MD_graphs/${name}_m2_${b2}
cp output_${dir2}/${b3} ./${name}_MD_graphs/${name}_m2_${b3}
cp output_${dir2}/${c1} ./${name}_MD_graphs/${name}_m2_${c1}
cp output_${dir2}/${c2} ./${name}_MD_graphs/${name}_m2_${c2}
cp output_${dir2}/${c3} ./${name}_MD_graphs/${name}_m2_${c3}
cp output_${dir2}/${d1} ./${name}_MD_graphs/${name}_m2_${d1}
cp output_${dir2}/${d2} ./${name}_MD_graphs/${name}_m2_${d2}
cp output_${dir2}/${d3} ./${name}_MD_graphs/${name}_m2_${d3}
cp output_${dir2}/${name}_m2_${e1} ./${name}_MD_graphs
cp output_${dir2}/${name}_m2_${e2} ./${name}_MD_graphs
cp output_${dir2}/${name}_m2_${e3} ./${name}_MD_graphs

echo "========== Moved all graphs into MD_graphs =========="


### Python analysis script for creating a graph report ###

sed -i "s/XpdbnameX/$name/g" graph_report.py
python3 ./graph_report.py

mkdir ${name}_MD_report 
mv *.png ${name}_MD_report

echo "Created graphical report for all MD simulations and replicas."
echo "========== Moved all PNG's into MD_report =========="
echo "Opening graphical reports:"

#xdg-open ${name}_MD_report/*_1_*
#xdg-open ${name}_MD_report/*_2_*
#xdg-open ${name}_MD_report/*_3_*
#xdg-open ${name}_MD_report/*_4_*
#xdg-open ${name}_MD_report/*_5_*


### Comparing original AlphaFold structures ###

echo "|------------------------------------------------------|"
echo "|= = = = = = = = = =  POST-ANALYSIS = = = = = = = = = =|" 
echo "|------------------------------------------------------|"

gmx rms -s ${dir0}/${name}_m0.pdb -f ${dir1}/${name}_m1.pdb -o ${name}_MD_graphs/${name}_rmsd_m0_m1_.xvg << EOF
0
0
EOF

gmx rms -s ${dir0}/${name}_m0.pdb -f ${dir2}/${name}_m2.pdb -o ${name}_MD_graphs/${name}_rmsd_m0_m2_.xvg << EOF
0
0
EOF

gmx rms -s ${dir1}/${name}_m1.pdb -f ${dir2}/${name}_m2.pdb -o ${name}_MD_graphs/${name}_rmsd_m1_m2_.xvg << EOF
0
0
EOF

mv $dir0/mdp_py_sh/rmsd_bars.py .
sed -i "s/XnameX/$name/g" rmsd_bars.py
python3 ./rmsd_bars.py

echo "|----------------------------------------------------|"
echo "|= = = = = = = = = =  MODEL RMSD  = = = = = = = = = =|" 
echo "|----------------------------------------------------|"
echo "> > > Expressed in nm < < <"

tail ${name}_models_rmsd.csv

mv ${name}_models_rmsd.csv ${name}_MD_graphs
mv ${name}_fig6_rmsd_bars.png ${name}_MD_report
#xdg-open ${name}_MD_report/${name}_rmsd_bars.png


### MMBPSA ###

echo "|------------------------------------------------|"
echo "|= = = = = = = = = =  MMPBSA  = = = = = = = = = =|" 
echo "|------------------------------------------------|"

# Model 0
mv  $dir0/mdp_py_sh/mmpbsa.in  $dir0
mv  $dir0/mdp_py_sh/mmpbsa_run.sh  $dir0
(cd $dir0 && ./mmpbsa_run.sh ${name}_m0)

# Model 1
mv  $dir1/mdp_py_sh/mmpbsa.in  $dir1
mv  $dir1/mdp_py_sh/mmpbsa_run.sh  $dir1
(cd $dir1 && ./mmpbsa_run.sh ${name}_m1)

# Model 2
mv  $dir2/mdp_py_sh/mmpbsa.in  $dir2
mv  $dir2/mdp_py_sh/mmpbsa_run.sh  $dir2
(cd $dir2 && ./mmpbsa_run.sh ${name}_m2)


mkdir MMPBSA
mv $dir0/MMPBSA_${name}_m0.dat MMPBSA
mv $dir1/MMPBSA_${name}_m1.dat MMPBSA
mv $dir2/MMPBSA_${name}_m2.dat MMPBSA
mv $dir0/MMPBSA_${name}_m0.csv MMPBSA
mv $dir1/MMPBSA_${name}_m1.csv MMPBSA
mv $dir2/MMPBSA_${name}_m2.csv MMPBSA

echo "========== Gathered all MMPBSA data =========="

echo "|-------------------------------------------|"
echo "| = = = = = = = =           = = = = = = = = |"
echo "| > > > > > > MASTER SCRIPT END < < < < < < |"
echo "| = = = = = = = =           = = = = = = = = |"
echo "|-------------------------------------------|"


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

