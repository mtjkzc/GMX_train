#!/bin/bash

###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################
#	  █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗
#	 ██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝
# 	 ███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗
#	 ██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║
#	 ██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║
#	 ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝                                                           
###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################

###  Analysis  script  of  MD  pipeline  ###

# 2024.06.13
# University of Zagreb, Facoulty of Science, Department of Chemistry
# Matej Kožić | mkozic@chem.pmf.hr

# master.sh runs MD script to run MD simulations
# anls_script.sh (this script)  follows after this script to analyze iz
# graph_report.py creates a graphical report after this script

# End of comment

###############################################################################
#}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}#
###############################################################################

# Run name
run_name=$1

# Defining stride
trjconv_skip=_trjconv_skip_


###  Finding  trajectory  file  in  current  directory  ###

directory="."
extension="_rawtraj.xtc"

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

name=${name/_rawtraj.xtc/}


###############################################################################
#[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][#
###############################################################################


###  Trajectory  processing  ###

echo "|---------------------------------------------------------------|"
echo "|= = = = = = = = = = CONVERTING TRAJECTORIES = = = = = = = = = =|"
echo "|---------------------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

gmx trjconv -f ${name}_rawtraj.xtc -s ${name}_md.tpr -center -ur compact -pbc nojump -skip $trjconv_skip -o ${name}.xtc << EOF
1
0
EOF


###  Analyses  ###

echo "|------------------------------------------------|"
echo "|= = = = = = = = = = ANALYSES = = = = = = = = = =|"
echo "|------------------------------------------------|"
echo "> > > Run name: "$run_name" < < <"

gmx rms -f ${name}.xtc -s ${name}_npt.gro -o ${name}_rmsd.xvg << EOF
1
1
EOF

gmx gyrate -f ${name}.xtc -s ${name}_npt.gro -o ${name}_gyr.xvg << EOF
1
EOF

gmx rmsf -res -f ${name}.xtc -s ${name}_npt.gro -o ${name}_rmsf.xvg << EOF
1
EOF

mkdir analysis_md
mv *.xvg analysis_md

echo 14 | gmx energy -f ${name}_md.edr -o analysis_md/md_temperature.xvg
echo 15 | gmx energy -f ${name}_md.edr -o analysis_md/md_pressure.xvg
echo 10 | gmx energy -f ${name}_md.edr -o analysis_md/md_potential.xvg

mkdir output_md
cp ${name}_md.tpr ${name}_md.edr ${name}_npt.gro ${name}.xtc output_md
cp analysis_md/* output_md

echo "|---------------------------------------------|"
echo "| = = = = = = = =             = = = = = = = = |"
echo "| > > > > > > ANALYSIS SCRIPT END < < < < < < |"
echo "| = = = = = = = =             = = = = = = = = |"
echo "|---------------------------------------------|"
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

