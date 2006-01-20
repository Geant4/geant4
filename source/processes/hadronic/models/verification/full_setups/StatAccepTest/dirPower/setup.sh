#
export DIR_INSTALLATIONS=/users/ribon/dirGrid/dirInstallations 
#
export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_BASE_DIR/lib
#
export PI_DIR=$DIR_INSTALLATIONS/dirPI 
export PATH=$PI_DIR/bin:${PATH} 
eval `aida-config --runtime sh` 
#
export GSL_DIR=$DIR_INSTALLATIONS/dirGSL 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib 
#
