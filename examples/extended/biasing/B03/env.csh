# This is an example script for setting the environment variables
# if you are on a redhat6.1 system at CERN.
# If you are running on a different platform
# you may have to install Anaphe5.0.1 and set the
# following variable ANAPHE_BASE appropriately.
# Please source this file before compiling and running
# the example.

#setenv ANAPHE_BASE /opt/Anaphe/5.0.1/specific/redhat72/gcc-2.95.2
setenv ANAPHE_BASE /afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2


setenv LD_LIBRARY_PATH ./:$G4WORKDIR/tmp/Linux-g++/exampleB03:$G4WORKDIR/lib/Linux-g++/:${LD_LIBRARY_PATH}
setenv PYTHONPATH $LD_LIBRARY_PATH
echo WARNING SETTING: G4LIB_BUILD_SHARED 1
setenv G4LIB_BUILD_SHARED 1
setenv PYTHON ${ANAPHE_BASE}/PublicDomainPackages/bin/python2.2
#setenv PATH /afs/cern.ch/sw/lhcxx/share/LHCXX/5.0.1/scripts/:$PATH
