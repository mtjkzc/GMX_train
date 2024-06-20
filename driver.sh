#!/bin/bash

###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################
#	 ██████╗ ██████╗ ██╗██╗   ██╗███████╗██████╗    ███████╗██╗  ██╗
#	 ██╔══██╗██╔══██╗██║██║   ██║██╔════╝██╔══██╗   ██╔════╝██║  ██║
#	 ██║  ██║██████╔╝██║██║   ██║█████╗  ██████╔╝   ███████╗███████║
#	 ██║  ██║██╔══██╗██║╚██╗ ██╔╝██╔══╝  ██╔══██╗   ╚════██║██╔══██║
#	 ██████╔╝██║  ██║██║ ╚████╔╝ ███████╗██║  ██║██╗███████║██║  ██║
#	 ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝  ╚══════╝╚═╝  ╚═╝╚═╝╚══════╝╚═╝  ╚═╝                                                               
###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################

# 2024.06.17 
# University of Zagreb, Facoulty of Science, Department of Chemistry
# Matej Kožić | mkozic@chem.pmf.hr 

# A script for sequentially running MD simulations starting from AlphaFold3
# zip file structures.

# Requirements to be in a dir same as this script:
#	- AlphaFold structures in their native zip file
#	- driver.sh (This script)
#	- monitor.sh
#	- MD_engine directory (cointains md_script and everythign else)

# End of comment

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

# Defining parameters
batch_run_name="Production_batch_Ugljan_5ns" #initiated jun 19 15:04
structures=( "o01" "o02" "o03" "o04" "o05" "o06" "o07" "o08" "o09" "o10" "o11" "o12" "o13" "o14" "o15" "o16" "o17" "o18" "o19" "o20" "o21" "o22" "o23" "o24" "o25" "o26" "o27" "o28" "o29" "o30" "o31" "o32" "o33" "o34" "o35" "o36" "o37" "o38" "o39" "o40" "o41" "o42" "o43" "o44" "o45" "o46" "o47" "o48" "o49" "o50" "o51" "o52" "o53" "o54" "o55" "o56" "o57" "o58" "o59" "o60" "o61" "o62" "o63" "o64" "o65" "o66" "o67" "o68" "o69" "o70" "o71" "o72" "o73" "o74" "o75" "o76" "o77" )


# Just copy from structures but between ' '
structure_string='( "o01" "o02" "o03" "o04" "o05" "o06" "o07" "o08" "o09" "o10" "o11" "o12" "o13" "o14" "o15" "o16" "o17" "o18" "o19" "o20" "o21" "o22" "o23" "o24" "o25" "o26" "o27" "o28" "o29" "o30" "o31" "o32" "o33" "o34" "o35" "o36" "o37" "o38" "o39" "o40" "o41" "o42" "o43" "o44" "o45" "o46" "o47" "o48" "o49" "o50" "o51" "o52" "o53" "o54" "o55" "o56" "o57" "o58" "o59" "o60" "o61" "o62" "o63" "o64" "o65" "o66" "o67" "o68" "o69" "o70" "o71" "o72" "o73" "o74" "o75" "o76" "o77" )'
	
###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

# Modifying monitor parameters
sed -i "s/XstructuresX/${structure_string}/g" monitor.sh
sed -i "s/Xbatch_run_nameX/$batch_run_name/g" monitor.sh

# Initializing logs
echo "= = = Initializing run $batch_run_name = = =" > ${batch_run_name}.log 	
echo "= = = MMPBSA report for $batch_run_name = = =" > mmpbsa_${batch_run_name}.log 
echo " " > mmpbsa_${batch_run_name}.log 

mkdir BATCH_RESULTS_$batch_run_name
mkdir MMPBSA_RESULTS_$batch_run_name
		
		
		
### Iterating over structures ###
		
for structure in "${structures[@]}"; do

	### MD simulation ###
	
	echo " " | tee -a ${batch_run_name}.log
	echo "|---------------------------------------------|" | tee -a ${batch_run_name}.log
	echo "|= = = = = = = = = =  ${structure}  = = = = = = = = = =|" | tee -a ${batch_run_name}.log
	echo "|---------------------------------------------|" | tee -a ${batch_run_name}.log
	echo " " | tee -a ${batch_run_name}.log

	mkdir $structure
	cp -r MD_engine/* $structure
	mv "fold_"${structure}".zip" $structure
	
	(cd $structure && ./master.sh $structure 2>&1 | tee -a ../${batch_run_name}.log )
	
	# Agregating results
	cp ${structure}/${structure}_af_MD_report/* BATCH_RESULTS_${batch_run_name}
	cp ${structure}/MMPBSA/* MMPBSA_RESULTS_$batch_run_name
	
	#######################################################################
	#/////////////////////////////////////////////////////////////////////#
	#######################################################################
	
	### MMPBSA ###
	
	echo "|---------------------------------------------|" >> mmpbsa_${batch_run_name}.log
	echo "|= = = = = = = = = =  ${structure}  = = = = = = = = = =|" >> mmpbsa_${batch_run_name}.log
	echo "|---------------------------------------------|" >> mmpbsa_${batch_run_name}.log
	
	grep "ΔTOTAL" "./MMPBSA_RESULTS_${batch_run_name}/MMPBSA_${structure}_af_m0.dat" >> mmpbsa_${batch_run_name}.log
	grep "ΔTOTAL" "./MMPBSA_RESULTS_${batch_run_name}/MMPBSA_${structure}_af_m1.dat" >> mmpbsa_${batch_run_name}.log
	grep "ΔTOTAL" "./MMPBSA_RESULTS_${batch_run_name}/MMPBSA_${structure}_af_m2.dat" >> mmpbsa_${batch_run_name}.log
	
	echo "= = = End of ${structure} = = =" >> mmpbsa_${batch_run_name}.log
	echo " " >> mmpbsa_${batch_run_name}.log
	
	echo "========== MMPBSA data appended ==========" | tee -a ${batch_run_name}.log
	echo " " | tee -a ${batch_run_name}.log
		
done


echo "|---------------------------------------|"
echo "| = = = = = = = =       = = = = = = = = |"
echo "| > > > > > DRIVER SCRIPT END < < < < < |"
echo "| = = = = = = = =       = = = = = = = = |"
echo "|---------------------------------------|"


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

