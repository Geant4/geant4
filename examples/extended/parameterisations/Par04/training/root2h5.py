"""
Converts a ROOT file to an HDF5 file, saving the shower energy in a 3D array.

Args:
    file_name: Path to the ROOT file with a name:
                output_<NAME>_angle<ANGLE>_<NUM>events_fullSim_<ID>.root
    output_dir: Path to the output directory, will create a file
                output_dir/<NAME>_Angle_<ANGLE>_<NUM>showers_<ID>.h5

Returns:
    None
"""

#!/bin/env python
import sys
import argparse
import numpy as np
import os
import uproot
import h5py

# Number of cells in the r, phi & z directions
NB_CELLS_R = 18
NB_CELLS_PHI = 50
NB_CELLS_Z = 45


def find_between(s, first, last):
    """
    Find a substring between two other substrings.

    Args:
        s (str): The string to search in.
        first (str): The first substring.
        last (str): The last substring.

    Returns:
        str: The substring found between 'first' and 'last'.
    """
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--outDir", type=str, default="")
    p.add_argument("--fName", type=str, default="")
    args = p.parse_args()
    return args


def main(argv):
    # Parse commandline arguments
    args = parse_args(argv)
    file_name = args.fName
    output_dir = args.outDir
    # The energy and angle of the particle are part
    # of the ROOT's file name so they can be extracted from the name
    # Get the energy value of the particle
    energy_particle = find_between(file_name, "output_", "_angle")
    # Get the angle value of the particule
    angle_particle = find_between(file_name, "_angle", "_")
    # Get the number of showers
    num_showers = find_between(file_name, "_", "events_")
    # Get Ids specific to the generated files
    file_id = find_between(file_name, "fullSim_", ".root")
    if os.stat(file_name).st_size > 0:
        h5_file = h5py.File(
            f"{output_dir}/{energy_particle}_Angle_{angle_particle}_{num_showers}showers_{file_id}.h5", "w"
        )
        # Read the Root file
        file = uproot.open(file_name)
        energy_particle = file["global"]["EnergyMC"].array()
        cell_r = file["virtualReadout"]["rhoCell"].array()
        cell_phi = file["virtualReadout"]["phiCell"].array()
        cell_energy = file["virtualReadout"]["EnergyCell"].array()
        cell_z = file["virtualReadout"]["zCell"].array()
        all_events = []
        # loop over events
        for event in range(len(energy_particle)):
            # Initialize a 3D array with shape nb_events, nb_cells in x,y,z
            data = np.zeros((NB_CELLS_R, NB_CELLS_PHI, NB_CELLS_Z))
            for ind in range(len(cell_r[event])):
                # This if statement is added to avoid having extra un-indexed cells
                if (
                    (cell_r[event][ind] < NB_CELLS_R)
                    and (cell_phi[event][ind] < NB_CELLS_PHI)
                    and (cell_z[event][ind] < NB_CELLS_Z)
                ):
                    data[cell_r[event][ind]][cell_phi[event][ind]][
                        cell_z[event][ind]
                    ] = cell_energy[event][ind]
            all_events.append(data)
        # Save dataset
        h5_file.create_dataset(
            f"{energy_particle}",
            data=all_events,
            compression="gzip",
            compression_opts=9,
        )
    h5_file.close()


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
