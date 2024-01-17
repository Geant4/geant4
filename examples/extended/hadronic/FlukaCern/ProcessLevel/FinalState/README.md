# Description

This example allows the simulation of hadron-nucleus inelastic nuclear interactions,    
and the study of the resulting final states.   
   
It is an adaptation of `Hadr09` example.  
It offers all `Hadr09` features, and adds the possibility 
of accessing hadron-nucleus inelastic nuclear interaction final states FROM `FLUKA`.  
  
With respect to the `Hadr09` example, the program also adds 
the possibility of PLOTTING the final state: 
secondaries energy spectra, and residual nuclei distributions.     
All plots (created via the G4 analysis manager) can be dumped 
to any of the usually supported formats (e.g. ROOT format), 
as well as in a Flair-compatible format.   
Regarding the extension of `G4H1` to insure `Flair` compatibility, 
see `geant4/examples/extended/hadronic/FlukaCern/utils`.
  
The class `HadronicGenerator` is the "generator".   
The main hadronic models (`CernFLUKAHadronInelastic`, `FTFP`, `QGSP`, `BERT`, `BIC`, `IonBIC`, `INCL`) are available.  
See `include/HadronicGenerator.hh` for more detailed information.    
In the `HadronicGenerator`, one can activate/desactivate coalescence and heavy fragments evaporation for `CernFLUKAHadronInelastic`.    
  
The main, `HadNucIneEvents.cc`, shows an example of how to use the event generator.    
In `HadNucIneEvents.cc`, one can select the physics models,
as well as the projectile hadron, its energy, its direction, the target material, and the number of collisions.
    
Note that the Geant4 run manager is not used.    
    
Before you can access the `FLUKA` hadron-nucleus inelastic models in this example, 
you will need to install and setup `FLUKA` and its interface.    
See the compulsory "Dependencies" paragraph below.    
   
A version of the interface to `FLUKA` is directly located at `geant4/examples/extended/hadronic/FlukaCern/FlukaInterface`.  
Note that for consistency, all calls to the random engine rely on the G4 random engine (including the calls from within the downloaded `FLUKA` release; see the `FlukaInterface` `GNUmakefile` to see how this is handled).   


# FLUKA inelastic hadron-nucleus interactions    
    
Hadron-NUCLEON interaction models are based on resonance production and decay below a few GeV,    
and on the Dual Parton model above.    
Hadron-NUCLEUS interactions: the PEANUT package includes    
a detailed Generalised Intra-Nuclear Cascade (GINC) and a preequilibrium stage,    
followed by equilibrium processes: evaporation, fission, Fermi break-up, gamma deexcitation.  
```  
A. Ferrari and P. Sala, “The Physics of High Energy Reactions,” in Proc. Workshop on Nuclear Reaction Data and Nuclear Reactors Physics, Design and Safety, p. 424, World Scientiﬁc, 1998.     
A. Ferrari and P. Sala, “Nuclear reactions in Monte Carlo codes,” Radiat. Prot. Dosimetry, vol. 99, no. 1-4, pp. 29–38, 2002.  
```  
      
   
# Dependencies

### Environment
- **gcc** >= 7 (Linux) and **gcc** >= 9 (MacOS)   
In practice, a recent version is recommended, at least `gcc >=10`.    
```
gcc --version
```

- **CMake** >= 3.16...3.21
```
cmake3 --version
```

- **G4** >= 11.0.3 (Not tested on older G4 releases: might still work, but with no guarantee).  
IMPORTANT: YOU NEED TO SOURCE YOUR G4 ENVIRONMENT.     
It needs to be sourced in whichever terminal you want to build / run a G4 application with the `FLUKA` interface.     
```
source path_to_geant4/install/bin/geant4.sh
which geant4-config   # NB: Your geant4-config should support the modern CMake way of building G4.
```

- **Easy setup on lxplus** (lxplus7):   
All you need to do on lxplus, to setup an environment satisfying all the conditions above, is, for example:
```
source /cvmfs/sft.cern.ch/lcg/releases/gcc/10.1.0/x86_64-centos7/setup.sh
source /cvmfs/geant4.cern.ch/geant4/11.1/x86_64-centos7-gcc10-optdeb-MT/CMake-setup.sh
# NB: Your geant4.sh is at: /cvmfs/geant4.cern.ch/geant4/11.1/x86_64-centos7-gcc10-optdeb-MT/bin/geant4.sh
```

### `FLUKA4`
Release: >= **4-3.2**     

Please install the latest `FLUKA` release.      
(1) You first need to register (and accept the licence when relevant): https://fluka.cern/download/registration   
(2) You can then download the `binary libraries` (or potentially the `source code` package, depending on your case):    
https://fluka.cern/download/latest-fluka-release.    
(3) Follow the `FLUKA` installation instructions: https://fluka.cern/documentation/installation    
In particular, for a Linux/MacOS install: https://fluka.cern/documentation/installation/fluka-linux-macos      
They will show you how to setup `FLUKA`.   
If (and only if) you went for the source code package option, you will need to build `fluka`, and, in addition, to do `make cpp_headers` at `path_to_fluka/src`.       
(4) Eventually, all you need are the headers `fluka_repo/include`, libraries `fluka_repo/lib`, and data `fluka_repo/data`. Check that they are not empty.    
Do not forget to add `/path_to_fluka/bin` to your `PATH`. Check with `which fluka`.  

### `FlukaInterface`
A version of the G4-FLUKA interface (`FLUKA` hadron-nucleus inelastic physics) 
is located at `geant4/examples/extended/hadronic/FlukaCern/FlukaInterface`.   
You will first need to build the interface to `FLUKA`, and create the environment scripts.   
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/FlukaInterface/
# Check with `which fluka` that fluka executable is added to your `PATH`.
$ source path_to_geant4/install/bin/geant4.sh
$ make interface
$ make env       # Creates `env_FLUKA.sh` and `env_FLUKA_G4_interface.sh`
```
IMPORTANT: `env_FLUKA_G4_interface.sh` needs to be sourced in whichever terminal 
you want to build / run a G4 application with the `FLUKA` interface.  


# Build this example
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/FinalState/
# Check with `which fluka` that fluka executable is added to your `PATH`.
$ source path_to_geant4/install/bin/geant4.sh
$ source ../../FlukaInterface/env_FLUKA_G4_interface.sh
$ mkdir build
$ cd build
$ cmake3 -DG4_USE_FLUKA=1 ..
$ make -j8
```


# Run this example
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/FinalState/
# Check with `which fluka` that fluka executable is added to your `PATH`.
$ source path_to_geant4/install/bin/geant4.sh
$ source ../../FlukaInterface/env_FLUKA_G4_interface.sh
$ ./build/HadNucIneEvents
```


# Study the final states
All plots are dumped at the end of the run in `all_secondaries.ext`.    
2 formats are supported: `ROOT` and `Flair`.


- You can use `ROOT`:
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/FinalState/
$ root all_secondaries.root
```

- Alternatively, the use of `Flair` is also supported.   
Please see http://flair.web.cern.ch/flair/download.html for `Flair` download. 
`Flair` tutorials are also available from that website.  
You can download the package corresponding to your distribution at the top of the page 
(no need for `geoviewer`, which is for geometry display). Then look at the requirements & installation instructions at the bottom of the page.  
If you face issues installing `Flair`, you can get support at: https://fluka-forum.web.cern.ch/c/installation/   
   
An example file, showing how to directly visualize the final states with Flair, is provided in this G4 example.   
By default, it directly provides comparison plots:  
`FTFP_BERT` versus `QGSP_BERT` versus `CernFLUKAHadronInelastic`, 7TeV proton on C.    
You can very easily adapt it to any study of interest (choice of physics models, choice of secondaries, etc).  
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/FinalState/
$ mkdir -p results/FLUKAHadronInelastic results/FTFP_BERT results/QGSP_BERT
$ 
$ # Choose physics case (modify HadNucIneEvents.cc), compile the G4 example, then run physics case:
$ cd results/FLUKAHadronInelastic
$ ../../build/HadNucIneEvents
$ # Etc for EACH physics case: FLUKAHadronInelastic, FTFP_BERT, QGSP_BERT.
$ 
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/FinalState/
$ ./update_final_state_flair_file.sh # Update `Det` indices in Flair file, to the ones observed in your simulation.
$ flair study_final_state.flair &
```
In the `Plot` tab, you can select the plot of interest in the left column, 
and then click `Plot` (top banner, yellow button).  
You can select a physics case by clicking on its name in the `Detectors` box (center). You can then decide to change its color, line width (`Options` box). You can decide to plot it or not, by selecting / unselecting `graph` in the `Show` box (in the center).    
IMPORTANT: You can select any secondary data (or residual nuclei data) which was created,
by chosing in the `Det` selection (button on the right).   
IMPORTANT: If a secondary does not appear in the `Det` drop-down menu, it means it is not part of the final state. In that case, you will want to unselect `graph`, so that no other secondary (see `Det`) is plotted.   
You can change the path of the data file by clicking on the folder button (button on the right).   
You can set the plots extrema, as well as select or unselect the log format, in the top right corner.  






