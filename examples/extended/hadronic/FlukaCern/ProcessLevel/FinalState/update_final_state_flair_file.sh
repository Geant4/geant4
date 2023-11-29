#!/usr/bin/env bash

# This script updates the `Det` fields in the Flair file, to the ones observed during the simulations.
# This is needed because the observed particles in the final state are DYNAMICALLY observed: 
# the number of particles changes from one simulation to the next, depending on the physics scenario, the number of events, etc.
# Hence, one needs to look what is the index of the `protons` histogram in the .hist file, etc.
# The use of this script is not needed in the XS G4 example case, because there, a FIXED number of XS are printed in the .hist files.

# Choose input files / number of comparison plots here
flair_file="study_final_state.flair"
num_comparison_plots=3
plot_block_max_length=100


# All plots
all_particles=$(cat $flair_file | grep "Plot:" | cut -d' ' -f2)

# Loop on all plots
for particle in ${all_particles[@]}; do

	echo "particle=$particle";
	
	# Line number of the "Plot: " block.
	plot_line_number=$(grep -n -m1 "Plot: $particle" $flair_file | cut -d':' -f1)
	echo "plot_line_number=$plot_line_number"
	
	# Find out if last plot block in the flair file.
	next_plot=$(grep -A$plot_block_max_length "Plot: $particle" $flair_file | grep -n -m2 "Plot: " | wc -l)
	if [ "$next_plot" -eq "2" ]; then
		has_next_plot=true
	else
		has_next_plot=false
	fi
	echo "has_next_plot=$has_next_plot";
	
	# Line number of the next "Plot: " block.
	if [ "$has_next_plot" = true ] ; then
		next_plot_line_number=$(grep -A$plot_block_max_length "Plot: $particle" $flair_file | grep -n -m2 "Plot: " | tail -n1 | cut -d':' -f1)
		# (next_plot_line_number - 1) is the extra number of lines taken by the plot block.
		next_plot_line_number=$(($plot_line_number + $next_plot_line_number - 1)) 
	else
		# Last line in the flair file.
		next_plot_line_number=$(cat $flair_file | wc -l)
	fi
	echo "next_plot_line_number=$next_plot_line_number"
	
	# Loop on all comparison plots.
	for ((i=0; i < $num_comparison_plots; i++)); do
		data_file=$(sed -n "$plot_line_number,$next_plot_line_number p" $flair_file | grep "file\.$i" | cut -d' ' -f2)
		echo "data_file=$data_file"
		# Data file not found: do nothing.
		if [ -z "$data_file" ]; then 
			echo "Warning: Tried to look for data file $data_file, which does not exist! Field is not updated for this file."
		# Found the data file.
		else
			# Get detector value in the data file.
			data_det=$(grep "# Detector:.* $particle" $data_file | cut -d' ' -f4)			
			data_det=$(($data_det - 1)) # An extra -1 because Flair detector indexing starts from 0.
			echo "data_det=$data_det"
	
			# Found detector value in data file.
			if [ ! -z "$data_det" ]; then
					
				# Look for detector field in the flair file.
				flair_det=$(sed -n "$plot_line_number,$next_plot_line_number p" $flair_file | grep "det\.$i")
				echo "flair_det=$flair_det"
			
				# Update value in detector field in flair file.
				if [ ! -z "$flair_det" ]; then
					sed -i "$plot_line_number,$next_plot_line_number s/det\.$i.*/det\.$i: $data_det/g" $flair_file
				# Flair file has no detector field in that plot: do not update field.
				else
					#sed -i "/file\.$i/i \ \tdet\.$i: $data_det" $flair_file
					echo "Warning: No det.$i defined for plot $particle, file $data_file: field was not updated."
				fi
				
			# Detector not found in data file: remove value in detector field in flair file.
			else
				sed -i "$plot_line_number,$next_plot_line_number s/det\.$i.*//" $flair_file
			fi
		fi
	done
	
done
