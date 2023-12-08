# Interface to `FLUKA` physics models

This is an interface to the `FLUKA` CERN hadron inelastic physics.  
It can be used by ANY G4 application (including the experimental frameworks, provided they are relying on G4 physics lists already).  
    
Note that this is a generic guide to be able to use the `FLUKA` interface from *any* G4 application.   
If you are trying to build / run the G4 examples included in `FlukaCern` with `FLUKA`,  
please directly follow the `README.md` files included in    
`ProcessLevel/CrossSection` and `ProcessLevel/FinalState`.  


# Description

The `FLUKA` hadron inelastic physics are accessible both through a `CrossSectionDataSet` and `HadronicInteraction`.  
The `CrossSectionDataSet` and `HadronicInteraction` can be used to define processes, themselves used to form a `FLUKAHadronInelasticPhysicsConstructor`.  
The `FLUKAHadronInelasticPhysicsConstructor` can in turn be used to form any `G4PhysicsList`. One example of `G4PhysicsList` is provided (initially `FTFP_BERT_HP`, but with hadron inelastic physics from `FLUKA` instead).  
Note that for consistency, all calls to the random engine rely on the G4 random engine (including the calls from within the downloaded `FLUKA` release; see the `FlukaInterface` `GNUmakefile` to see how this is handled).


# Structure

The project contains:  
- `fluka4_wrapper`: wrapper to `FLUKA`.
- `cpp_fortran_bridges`: fortran <-> C++ interfacing tools.
- `cpp_utils`: general project helpers (C++ string manipulation, etc).
- `fluka5`: interface to `FLUKA` inelastic hadron-nucleus interactions.
- `fluka_G4_physics_list`: `ModularPhysicsList`, `PhysicsConstructor`, `CrossSectionDataSet`, `HadronicInteraction` to make up a G4 physics list.
- `fluka_G4_bridges`: conversion tools to switch between `FLUKA` and G4 worlds (particles, random numbers, etc).


# Dependency

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
you want to compile / run a G4 application with the `FLUKA` interface.    


# Use as a FLUKA physics list (in ANY G4 application)
Reminder: this is a generic guide to be able to use the `FLUKA` interface from *any* G4 application.  
If you are trying to build / run the G4 examples included in `FlukaCern` with `FLUKA`, you can directly follow the `README.md` files included in those examples.  
Notably, the file in `cmake/Modules` is already included in those examples (the rest of the procedure is the same).  
  
```bash
$ cd G4_application
```

- Edit your G4 application's `CMakeLists.txt`:  
You need to add the `FLUKA` module dependency. You need to add the following (for example, after the `project` line):
```bash
#----------------------------------------------------------------------------
# FindFLUKAInterface.cmake is located at your_path_to_geant4/cmake/Modules/FindFLUKAInterface.cmake
# Check that FindFLUKAInterface.cmake can be found from $CMAKE_MODULE_PATH
message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
# Otherwise, you can always prepend it to the cmake module search path with:
# set(CMAKE_MODULE_PATH my_path_to_find_fluka ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------------
# Check whether FLUKA should be used or not
set(G4_USE_FLUKA OFF CACHE BOOL "Using FLUKA")
if(G4_USE_FLUKA)
  message(STATUS "G4_USE_FLUKA=ON : Using FLUKA interface for building ${PROJECT_SOURCE_DIR}")
  add_definitions(-DG4_USE_FLUKA)
  find_package(FLUKAInterface REQUIRED)
  if(FLUKAInterface_FOUND)
    message(STATUS "FLUKA cmake module was found : ${CMAKE_MODULE_PATH}")
  else()
    message(FATAL_ERROR "FLUKA cmake module was NOT found! Please add one.")
  endif()
else()
  message(STATUS "G4_USE_FLUKA=OFF : NOT using FLUKA interface for building ${PROJECT_SOURCE_DIR}. \n \
  If ever you want to use the FLUKA interface, please repeat cmake command with -DG4_USE_FLUKA=1")
endif()
```
Additionally, add:  
`${FLUKAInterface_INCLUDE_DIR}` to the `include_directories` line.  
`${FLUKAInterface_LIBRARIES}` to the `target_link_libraries` line.  

- You can now compile and run your G4 application:
```bash
# Check with `which fluka` that fluka executable is added to your `PATH`.
$ source path_to_geant4/install/bin/geant4.sh
$ source path_to_FlukaInterface/env_FLUKA_G4_interface.sh             # Carefully follow the setup instructions if anything fails.
$ cd build
$ cmake3 -DGeant4_DIR=your_path_to_geant4 -DG4_USE_FLUKA=1 ../
$ make -j8
$ ./G4_executable
```

