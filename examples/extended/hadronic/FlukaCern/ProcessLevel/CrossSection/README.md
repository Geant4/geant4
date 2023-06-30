# Description

This example allows the study of G4 cross-sections, 
and in addition, of the `FLUKA` hadron-nucleus reaction cross sections.

The user can printout any particle-material XS.    
The XS are exactly the ones defined in any G4 physicsList chosen by the user,
or from `FLUKA` (hadron-nucleus inelastic case).
  
In the input file, the user can set:
- projectile.
- target material (element, compound or even mixture).
- plotting options.

All plots (created via the G4 analysis manager) can be dumped 
to any of the usually supported formats (e.g. `ROOT` format), 
but also in a Flair-compatible format.   
Regarding the extension of `G4H1` to insure `Flair` compatibility, 
see `geant4/examples/extended/hadronic/FlukaCern/utils`.
     
Note that the Geant4 run SERIAL manager is used, since `FLUKA` is single-threaded  
(in actual `FLUKA` runs, parallelism is achieved via a multi-processing approach).    
   
Before you can access the `FLUKA` hadron-nucleus inelastic model in this example, 
you will need to install and setup `FLUKA` and its interface.    
See "Dependencies" paragraph below.    

A version of the interface to `FLUKA` is directly located at `geant4/examples/extended/hadronic/FlukaCern/FlukaInterface`.      
   

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
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/CrossSection/
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
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/CrossSection/
# Check with `which fluka` that fluka executable is added to your `PATH`.
$ source path_to_geant4/install/bin/geant4.sh
$ source ../../FlukaInterface/env_FLUKA_G4_interface.sh
# Edit all_XS.in (choice of particle, material, energy range...)
# XS used by FTFP_BERT_HP physics list:
$ ./build/HadronNucleusXS all_XS.in FTFP_BERT_HP
# XS used by G4_HP_CernFLUKAHadronInelastic_PhysicsList physics list:
$ ./build/HadronNucleusXS all_XS.in G4_HP_CFLUKAHI
```


# Study the cross-sections
All plots are dumped at the end of the run in `all_XS.ext`.    
2 formats are supported: `ROOT` and `Flair`.


- You can use `ROOT`:
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/CrossSection/
$ root all_XS.root
```

- Alternatively, the use of `Flair` is also supported.   
Please see http://flair.web.cern.ch/flair/download.html for `Flair` download. 
`Flair` tutorials are also available from that website.  
You can download the package corresponding to your distribution at the top of the page 
(no need for `geoviewer`, which is for geometry display). Then look at the requirements & installation instructions at the bottom of the page.  
If you face issues installing `Flair`, you can get support at: https://fluka-forum.web.cern.ch/c/installation/  
  
An example file, showing how to directly visualize all XS with `Flair`, 
is provided in this G4 example.   
By default, it directly provides comparison plots:  
`FTFP_BERT_HP` versus `G4_HP_CernFLUKAHadronInelastic`.    
One can very easily adapt it (directly from `Flair` GUI) to any study of interest.
```bash
$ cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/CrossSection/
$ mkdir -p results/FTFP_BERT_HP results/G4_HP_CernFLUKAHadronInelastic
$ # First run the example for FTFP_BERT, G4_HP_CernFLUKAHadronInelastic, and move the results!
$ ./build/HadronNucleusXS all_XS.in FTFP_BERT_HP
$ mv all_XS.hist results/FTFP_BERT_HP
$ ./build/HadronNucleusXS all_XS.in G4_HP_CFLUKAHI
$ mv all_XS.hist results/G4_HP_CernFLUKAHadronInelastic
$ 
$ flair study_all_XS.flair &
```
In the `Plot` tab, you can select the plot of interest in the left column, 
and then click `Plot` (top banner, yellow button).    
You can select a physics list by clicking on its name in the `Detectors` box (center). You can then decide to change its color, line width (`Options` box). You can decide to plot it or not, by selecting / unselecting `graph` in the `Show` box (in the center).  
You can change the path of the data file by clicking on the folder button (button on the right).   
You can set the plots extrema, as well as select or unselect the log format, in the top right corner.  



