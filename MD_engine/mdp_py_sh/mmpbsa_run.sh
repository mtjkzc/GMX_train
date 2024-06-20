#!/bin/bash

replica_name=$1

source /opt/amber/amber.sh
source /opt/gromacs/bin/GMXRC.bash
source /opt/mmpbsa/bin/activate

gmx_MMPBSA -O -i mmpbsa.in -cs ${replica_name}_md.tpr -ct ${replica_name}.xtc -ci ${replica_name}_CHAINS.ndx -cg 2 3 -cp  ${replica_name}.top -o MMPBSA_${replica_name}.dat -eo MMPBSA_${replica_name}.csv -nogui
