#
# --- Setup the Geant4 kernel environment ---
# ***LOOKHERE*** IT IS ASSUMED THAT YOU HAVE ALREADY YOUR 
#                COMPILER g++ 3.23 AND THE GEANT4 ENVIRONMENTAL 
#                VARIABLES PROPERLY SET.
###source ~/ExtraSpace/geant4/env.csh-3.23-nodeb
#
# --- Setup your Geant4 application environment ---
unsetenv G4TMP
unsetenv G4BIN
#
#***LOOKHERE*** SET HERE YOUR WORKSPACE: BY DEFAULT, IT IS
#               ASSUMED TO BE THE CURRENT DIRECTORY.
setenv G4WORKDIR $PWD
#
setenv PATH ${PATH}:${G4WORKDIR}/bin/${G4SYSTEM}
setenv G4ANALYSIS_USE 1
#
# --- Aida / PI ---
setenv PI_DIR /afs/cern.ch/sw/lcg/app/releases/PI/PI_1_2_3/rh73_gcc323
setenv PATH ${PI_DIR}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PI_DIR}/lib
eval `aida-config --runtime csh`
#
# --- GSL : this is needed only for  dirStat/pvalue.cpp ---
setenv GSL_DIR /afs/cern.ch/sw/lcg/external/GSL/1.4/rh73_gcc323
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GSL_DIR}/lib
#

