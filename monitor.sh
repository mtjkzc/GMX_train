#!/bin/bash

###############################################################################                                                  
# START START START START START START START START START START START START STA #
###############################################################################
#	 ███╗   ███╗ ██████╗ ███╗   ██╗██╗████████╗ ██████╗ ██████╗ 
#	 ████╗ ████║██╔═══██╗████╗  ██║██║╚══██╔══╝██╔═══██╗██╔══██╗
#	 ██╔████╔██║██║   ██║██╔██╗ ██║██║   ██║   ██║   ██║██████╔╝
#	 ██║╚██╔╝██║██║   ██║██║╚██╗██║██║   ██║   ██║   ██║██╔══██╗
#	 ██║ ╚═╝ ██║╚██████╔╝██║ ╚████║██║   ██║   ╚██████╔╝██║  ██║
#	 ╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝
###############################################################################                                                  
# START START START START START START START START START START START START STA #
###############################################################################                                                         

# 2024.06.18 
# University of Zagreb, Facoulty of Science, Department of Chemistry
# Matej Kožić | mkozic@chem.pmf.hr 

# Monitors the progress of batch simulations that are run by driver.sh.
# The script iterates over inverted array of structures that driver.sh runs.
# It looks for the latest log to append to progress_log and tail it on display.
# The script has to be run in a separate terminal, while driver.sh is running.
# It has to be ran manually.
# To run it continuous run with a flag "-c" ( ./monitor.sh -c )

# End of comment

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

# Initial parameters (driver.sh assigns this upon initialization)
batch_run_name="Xbatch_run_nameX"
structure_list=XstructuresX

update_delay=1 # seconds
clear_after=20 # seconds
continuous=$1

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

### Reversing the order of structures ###
 
reversed_s_list=()

for ((i=${#structure_list[@]}-1; i>=0; i--)); do
	reversed_s_list+=("${structure_list[i]}")
done


# Initializing global log
#echo "= = = Tracking progress for ${batch_run_name} = = =" > ./progress.log

# Defining file structure
structures=$reversed_s_list
replica_directories=( "replica_2" "replica_1" "replica_0" )
logs=( "md_output.log" "npt_output.log" "nvt_output.log" "min_output.log" )


### Engine loop ###

counter=1
clear_counter=0
while true; do
	echo " "
    if [ $counter -eq 1 ]; then
        echo " = = = ( < ( < ( < ( < ( < ( < ( = = ="
        counter=0
        ((clear_counter++))
    else
        echo " = = = > ) > ) > ) > ) > ) > ) > = = = "
        counter=1
        ((clear_counter++))
    fi
    if [ ${clear_counter} = ${clear_after} ]; then
    	clear
    fi
	for structure in "${reversed_s_list[@]}"; do
		if [ -e $structure ]; then
			for rdir in "${replica_directories[@]}"; do
				if [ -e ${structure}/${rdir} ]; then
					for log in "${logs[@]}"; do
						if [ -e ${structure}/${rdir}/${log} ]; then
							
							echo " > > > Run name: $batch_run_name"
							echo " > > > Queue: ${structure_list[@]}"
							echo " > > > ${structure} | ${rdir} | ${log} < < <" #| tee -a ./progress.log
							echo " " #| tee -a ./progress.log
							tail ${structure}/${rdir}/${log} #>> ./progress.log
							echo " " #| tee -a ./progress.log
							
							if [ "${continuous}" = "-c" ]; then
								break 3
							else
								break 4
							fi
						fi
					done
				fi
			done
		fi
	done
	sleep $update_delay
done


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

