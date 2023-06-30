# stim_pixe_tomography advanced example

The stim_pixe_tomography advanced example is developed to simulate three dimensional STIM or
PIXE tomography experiments. The simulation results are written in a binary file and can be easily accessed using the
provided scripts.

 Publications:

    [1] Li Z, Incerti S, Beasley D, Shen H, Wang S, Seznec H, et al. Accuracy of three-dimensional proton imaging of an
    inertial confinement fusion target assessed by Geant4 simulation. Nucl Instrum Methods Phys Res B. 2023;
    536:38-44. https://doi.org/10.1016/j.nimb.2022.12.026.

    [2] Michelet C, Li Z, Jalenques H, Incerti S, Barberet P, Devs G, et al. A Geant4 simulation of X-ray emission
    for three-dimensional proton imaging of microscopic samples. Phys Med. 2022;94:85-93. https://doi.org/10.1016/j.ejmp.2021.12.002.

    [3] Michelet C, Li Z, Yang W, Incerti S, Desbarats P, Giovannelli J-F, et al. A Geant4 simulation for
    three-dimensional proton imaging of microscopic samples. Phys Med. 2019;65:172-80. https://doi.org/10.1016/j.ejmp.2019.08.022.

 Contact:

    michelet@lp2ib.in2p3.fr      (Claire Michelet)

    zhuxin.li@outlook.com        (Zhuxin Li)


 More information and a detailed UserGuide are available:
http://geant4.in2p3.fr  (Documentation section)

## 1 - GEOMETRY DEFINITION

  Three phantoms are available, users can build up new phantoms or choose the following
  three phantoms by setting the "phantom_type":

    1) A simple cube (see publication [2-3]), phantom_type = 1

    The absorber is a box made of a given material.

    2) Upper part of Caenorhabditis elegans (C.elegans) worm (see publication [2-3]) , phantom_type = 2

    C.elegans phantom is composed of  6 ellipsoids. The size and shape of ellipsoids are based on the
    nanotoxicology studies carried-out at LP2I Bordeaux laboratory .

    3) Inertial confinement fusion (ICF) target (see publication [1]), phantom_type = 3

    ICF target is sphere shell, made of Ge-doped glow discharge polymer (GDP)


 ##2 - PHYSICS LIST

  Physics lists are based on modular design. Several modules are instantiated:

    1) Transportation
    2) EM physics
    3) Decay physics
    4) Hadron physics, optional

 EM physics builders can be local or from G4 kernel physics_lists subdirectory.

    - "emlivermore"             default low-energy EM physics using Livermore data
    - "local"                   local physics builders, options are explicit in PhysListEmStandard
    - "emstandard_opt0"         recommended standard EM physics for LHC
    - "emstandard_opt1"         best CPU performance standard physics for LHC
    - "emstandard_opt2"         similar fast simulation
    - "emstandard_opt3"         best standard EM options - analog to "local" above
    - "emstandard_opt4"         best current advanced EM options standard + lowenergy
    - "emstandardWVI"           standard EM physics and WentzelVI multiple scattering
    - "emstandardSS"            standard EM physics and single scattering model
    - "emstandardGS"            standard EM physics and Goudsmit-Saunderson multiple scatt.
    - "empenelope"              low-energy EM physics implementing Penelope models
    - "emlowenergy"             low-energy EM physics implementing experimental

 Decay and StepMax processes are added to each list.

 Optional components can be added:

    - "elastic"       elastic scattering of hadrons
    - "binary"        QBBC configuration of hadron inelastic models
    - "binary_ion"    Binary ion inelastic models

 Physics lists and options can be (re)set with UI commands.


  ##3 - HOW TO RUN

 To run a PIXE tomography simulation in 'batch' mode using a pixe3d.mac file:

    ./stim_pixe_tomography -p pixe3d.mac

 or if you want to specify the number of threads:

    ./stim_pixe_tomography -p pixe3d.mac N

    N is the number of threads

  An example of pixe3d.mac is provided.
  It is designed for the PIXE-T simulation of the cube phantom of 40 um.
  It is defined for 10 projections  1 slice  20 pixels. 1000000 protons are used for each beam.


 To run a STIM tomography simulation:

      ./stim_pixe_tomography -s pixe3d_stim.mac

  or if you want to specify the number of threads:

    ./stim_pixe_tomography -s pixe3d_stim.mac N

    N is the number of threads

 An example of pixe3d_stim.mac (arbitrarily name, you may rename it pixe3d.mac if you wish) is provided.
 It is designed for the STIM-T simulation of the cube phantom of 40 um.
 It is defined for 10 projections  1 slice  20 pixels. 100 protons are used for each beam.

 ##4 - VISUALISATION
To visualize the phantoms, run:

    ./stim_pixe_tomography

 ##5 - OUTPUT FILES

 If a PIXE tomography simulation is made, two files are going to be generated:

    1) GammaAtCreation.dat, which keeps the info of secondary photons at creation
    2) GammaAtExit.dat, which keeps the info of secondary photons at exit of the phantom

 If a STIM tomography simulation is made, ProtonAtExit.dat is generated, in which the info of primary protons is kept


 ##6 - LIST OF MACROS AND SCRIPTS

Once you build the example, the following macros and script will be copied to your build directory:

       pixe3d.mac: an example macro to run a PIXE-T simulation for cube of 40 um
       pixe3d_stim.mac: an example macro to run a STIM-T simulation for cube of 40 um
       pixe3d_initial.mac: it contains the information of physics processes
       init_vis.mac and vis.mac: for the visualization
       GPSPointLoop.C: it generates a macro file to run the simulation by reading pixe3d_initial.mac

In the **Scripts** folder, you will find other scripts for different uses.

***To obtain the reconstruction data:***

       BinToStd_ProtonAtExit.C: it reads the STIM-T simulation results and generates the data file for STIM-T reconstruction using selection with particle momentum.
       BinToStd_GammaAtCreation.C: it reads the PIXE-T simulation results for X-rays at creation and generates the data file for PIXE-T reconstruction using selection with particle momentum.
       BinToStd_GammaAtExit.C: it reads the PIXE-T simulation results for X-rays at exit and generates the data file for PIXE-T reconstruction using selection with particle momentum.
       BinToStd_proton_position.C: it reads the STIM-T simulation results and generates the data file for STIM-T reconstruction using selection with particle position and momentum
       BinToStd_gamma_position.C: it reads the PIXE-T simulation results for X-rays and generates the data file for PIXE-T reconstruction using selection with particle position and momentum

 ***To locate the interruption if an interruption of simulation occurs:***

       LocateInterruption_ProtonAtExit.C: in case of interruption, it locates the projection position of interruption for STIM-T simulation.
       LocateInterruption_GammaAtExit.C: in case of interruption, it locates the projection position of interruption for PIXE-T simulation.
***To obtain the reconstruction data in case of an interruption of simulation:***

       Concatenate_BinToStd_ProtonAtExit.C: in case of one interruption, it reads STIM-T simulation results and generates the data file for STIM-T reconstruction.
       Concatenate_BinToStd_GammaAtCreation.C: in case of one interruption, it reads PIXE-T simulation results for X-rays at creation and generates the data file for PIXE-T reconstruction.
       Concatenate_BinToStd_GammaAtExit.C: in case of one interruption, it reads PIXE-T simulation results for X-rays at exit and generates the data file for PIXE-T reconstruction.
***To visualize the spectrum:***

       Spectrum_proton.C: it visualizes the spectrum of protons and plots a histogram by reading simulation result ProtonAtExit.dat.
       Spectrum_gamma.C: it visualizes the spectrum of X-rays and plots a histogram by reading simulation result GammaAtCreation.dat or GammaAtExit.dat.
       TomoSpectrum_HIST_proton.C: it visualizes the spectrum of protons and plots a histogram by reading StimEvent data. It also writes the spectrum data in a txt file.
       TomoSpectrum.C: it visualizes the spectrum of X-rays and plots a graph by reading PixeEvent data. It also writes the spectrum data in a txt file.
       TomoSpectrum_HIST.C: it visualizes the spectrum of X-rays and plots a histogram by reading PixeEvent data. It also writes the spectrum data in a txt file.
***Scripts for specific use:***

       Extract_Projection.C: it extracts 50 projections from a PixeEvent data file for tomographic reconstruction, which contains 100 projections. In fact, it extracts the projection 0, 2, 4, 6, 898 from projections 0-99. It eventually generates a new file with new index number of projections 0-49.
       Check_PixeEventFile.C: it checks if the index of projections of a PixeEvent data file for tomographic reconstruction is correct. For example, if the user extract 50 projections from a data file composed 100 projections, it is necessary to make sure in the new data file, the index of projection starts from 0 and ends at 49.
       Extract_Slice.C: it extracts a certain number of slice(s) from a PixeEvent data file for tomographic reconstruction. Users need to specify the first and the last slice to be extracted. Note that when writing a new data file, the index of slices will be initiated from 0.
       Concatenate_BinToStd_GammaAtCreation_fabricate.C: if users make a PIXE-T simulation on a symmetrical object with only one projection, this script can be used to fabricate the other 99 projection data for X-rays at creation with same energy.
       Concatenate_BinToStd_GammaAtExit_fabricate.C: if users make a PIXE-T simulation on a symmetrical object with only one projection, this script can be used to fabricate the other 99 projection data for X-rays at exit with same energy
***Scripts to generate voxelized phantoms:***

In order to compare the reconstructed tomographic images with original
phantoms, it may be necessary to use a voxelized phantom.

       generate_voxelized_sphere_phantom.py: it generates a voxelized phantom of an inertial confinement fusion target.
       generate_voxelized_worm_phantom.py: it generates a voxelized phantom of the upper part of C. elegans.
More information can be found in the UserGuide.