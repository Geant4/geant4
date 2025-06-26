"""
Converts an EDM4HEP ROOT file to an HDF5 file, saving the shower energy in a 3D array and shower data in a 2D array.
Units: Energy values are stored in MeV and angles are stored in radians.

"""

#!/bin/env python
import sys
import argparse
import numpy as np
import os
import uproot
import h5py


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--outputDir", '-o', type=str, default="./", help="Path to the output directory")
    p.add_argument("--inputFile", '-i', type=str, required=True, help="Name of the EDM4hep file to translate")
    p.add_argument("--numR", type=int, default=18, help="Number of cells in R")
    p.add_argument("--numPhi", type=int, default=50, help="Number of cells in phi")
    p.add_argument("--numZ", type=int, default=45, help="Number of cells in z")
    p.add_argument("--samplingFraction", type=float, default=1., help="Sampling fraction to use to rescale cell energy. Defined as f=active/(active+absorber)")
    args = p.parse_args()
    return args


def main(argv):
    # Parse commandline arguments
    args = parse_args(argv)
    input_file = args.inputFile
    output_dir = args.outputDir

    # Number of cells in the r, phi & z directions
    num_cells_R = args.numR
    num_cells_phi = args.numPhi
    num_cells_z = args.numZ

    # Sampling fraction that rescales energy of each cell
    sampling_fraction = args.samplingFraction

    if os.stat(input_file).st_size > 0:
        h5_file = h5py.File(
            f"{output_dir}/{os.path.splitext(os.path.basename(input_file))[0]}.h5", "w"
        )
        print(f"Creating output file {output_dir}/{os.path.splitext(os.path.basename(input_file))[0]}.h5")
        # Read Root file
        file = uproot.open(input_file)
        energy_particle = file["global"]["EnergyMC"].array()
        # For future once theta,phi are implemented in Par04 event/run action
        #phi_particle = file["global"]["PhiMC"].array()
        #theta_particle = file["global"]["ThetaMC"].array()
        cell_r = file["virtualReadout"]["rhoCell"].array()
        cell_phi = file["virtualReadout"]["phiCell"].array()
        cell_energy = file["virtualReadout"]["EnergyCell"].array()
        cell_z = file["virtualReadout"]["zCell"].array()
        all_events = []
        num_showers = len(energy_particle)
        # loop over events
        for event in range(num_showers):
            # Initialize a 3D array with shape nb_events, nb_cells in x,y,z (rho,phi,z)
            shower = np.zeros((num_cells_R, num_cells_phi, num_cells_z))
            for cell in range(len(cell_r[event])):
                # This if statement is added to avoid having cells outside of desired cylinder size
                if (
                    (cell_r[event][cell] < num_cells_R)
                    and (cell_phi[event][cell] < num_cells_phi)
                    and (cell_z[event][cell] < num_cells_z)
                ):
                    shower[cell_r[event][cell]][cell_phi[event][cell]][
                        cell_z[event][cell]
                    ] = cell_energy[event][cell]
            all_events.append(shower)
        # Save dataset
        print(f"Creating datasets with shape {np.shape(energy_particle)} and {np.shape(all_events)} ")
        h5_file.create_dataset("incident_energy", data=energy_particle, compression="gzip", compression_opts=9,)
        # For future once theta,phi are implemented in Par04 event/run action
        #h5_file.create_dataset("incident_phi", data=phi_particle, compression="gzip", compression_opts=9,)
        #h5_file.create_dataset("incident_theta", data=theta_particle, compression="gzip", compression_opts=9,)
        h5_file.create_dataset("showers", data=all_events, compression="gzip", compression_opts=9,)
    h5_file.close()


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
