setenv LD_LIBRARY_PATH ./:$G4WORKDIR/tmp/Linux-g++/exampleB03:$G4WORKDIR/lib/Linux-g++/:${LD_LIBRARY_PATH}
setenv PYTHONPATH $LD_LIBRARY_PATH
echo WARNING SETTING: G4LIB_BUILD_SHARED 1
setenv G4LIB_BUILD_SHARED 1
setenv PYTHON /afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/PublicDomainPackages/bin/python2.2
