# This is an example script for setting the environment variables
# if you are on a redhat6.1 system at CERN.
# Please source this file before compiling and running
# the example.

#if using afs
setenv PYTHON_DIR /afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/PublicDomainPackages/2.0.0

# if Anaphe is installed in /opt
#setenv PYTHON_DIR /opt/Anaphe/5.0.1/specific/redhat72/gcc-2.95.2/PublicDomainPackages/2.0.0


setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}:${G4WORKDIR}/tmp/${G4SYSTEM}/exampleB03
setenv PYTHONPATH $LD_LIBRARY_PATH.

echo WARNING SETTING: G4LIB_BUILD_SHARED 1
setenv G4LIB_BUILD_SHARED 1

