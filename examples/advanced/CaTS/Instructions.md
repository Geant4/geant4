# Prerequisites
Opticks requires Geant4 (10.7.p02), nvidia cuda (11.3)  and nvidia Optix (6.5) among other libraries. CaTS in addition will require ROOT. If all these libraries and development headers are available on your machine skip directly to  (**Building opticks vs. existing libraries**). On a 'blank' computing system it makes sense to build CLHEP, then Geant4 and finally ROOT assuring that all the necessary development libraries and headers are installed.   

# Building CLHEP
The current version of Geant4 10.07.p02 is build on clhep 2.4.4.0. 
CLHEP can be found at:
https://proj-clhep.web.cern.ch/proj-clhep/clhep23.html

to build it from scratch using cmake (used cmake version 3.20.5) 

    cd to the directory where you want to build clhep
    wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep-2.4.4.0.tgz
    tar xzvf clhep-2.4.4.0.tgz
    cd 2.4.4.0/
    mkdir CLHEP-build
    cd  CLHEP-build
    cmake -DCMAKE_INSTALL_PREFIX=../CLHEP-install DCLHEP_BUILD_CXXSTD=-std=c++11 ../CLHEP
    make -j 8
    make install

**Note** the default install directory is /usr/local but one needs root privileges to install it there:

    cd to the directory where you want to build clhep
    wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep-2.4.4.0.tgz
    tar xzvf clhep-2.4.4.0.tgz
    cd 2.4.4.0/
    mkdir CLHEP-build
    cd  CLHEP-build
    cmake -DCLHEP_BUILD_CXXSTD=-std=c++11 ../CLHEP
    make -j 8
    sudo make install

# Building Geant4

Geant4 versions are available at:
https://geant4.web.cern.ch/support/download


    cd to the directory where you want to install Geant4
    wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p02.tar.gz
    mkdir geant4.10.07.p02-build
    cd  geant4.10.07.p02-build
    cmake -DCMAKE_INSTALL_PREFIX=../geant4.10.07.p02-install -DGEANT4_BUILD_VERBOSE_CODE=OFF -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_SYSTEM_EXPAT=ON -DGEANT4_USE_SYSTEM_ZLIB=ON  -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_USE_OPENGL_X11=ON ../geant4.10.07.p02
    make -j 8
    make install
    . ../geant4.10.07.p02-install/bin/geant4.sh


check the output for any error, install any development packages that might be necessary. 




# Building ROOT 
Instructions how to build ROOT from ssource can be found at:
https://root.cern/install/build_from_source/

    # cd  to the diretory where you want to install root
    git clone --branch latest-stable https://github.com/root-project/root.git root_src
    mkdir root-build
    cd root-build
    cmake -DCMAKE_INSTALL_PREFIX=../root-install  ../root
    # speed up the make process
    new=" -j$(($(grep -c ^processor /proc/cpuinfo) - 1))" 
    case ":${MAKEFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) MAKEFLAGS="$MAKEFLAGS:$new"  ;;
    esac
    cmake -DCMAKE_INSTALL_PREFIX=../root-install  ../root_src/
    cmake --build . --target install

check the output for any error, install any development packages that might be necessary. 


    . ../root-install/bin/thisroot.sh
    root
    

# Installing CUDA

cuda (11.3) is available at the NVIDIA web site just follow the instruction depending on the system you are using. 

https://developer.nvidia.com/cuda-downloads

**Note** this will also install the corresponding NVIDIA graphics driver you might have to reboot.  

Another way to obtain cuda is to install the NVIDIA hpc-sdk kit. Which can be found here. 
https://developer.nvidia.com/nvidia-hpc-sdk-downloads

The NVIDIA hpc-sdk kit provides an interesting set of tools e.g. nvc++ which allows offloading of parallel algorithms to NVIDIA GPUs.

A good way to check that things are working properly is to build the cuda samples and execute them

    # cd to the directory where you want to build the cuda samples. E.g. the commands deviceQueryDrv and deviceQuery provide useful information. 
    mkdir cuda-test
    cd cuda-test
    cp -r /usr/local/cuda-11.3/samples .
    cd  samples/
    which nvcc
    make 
    bin/x86_64/linux/release/deviceQuery
    bin/x86_64/linux/release/deviceQueryDrv 

2 Tools For Monitoring Nvidia GPUs On Linux can be found here:
https://www.linuxuprising.com/2019/06/2-tools-for-monitoring-nvidia-gpus-on.html

# Installing Optix (6.5)

https://developer.nvidia.com/designworks/optix/download

Optix comes with precompiled samples and one might want to try them:

    # cd to the Optix installation directory
    cd SDK-precompiled-samples
    export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
    # execute e.g.:
    ./optixMDLSphere
    ./optixSphere
    ./optixTutorial
    # etc.





# Building opticks vs. existing libraries

This are instructions how to build opticks making use of preinstalled libraries available on the system. These libraries include CLHEP, xerces-c, boost and  Geant4.
For geant 4 we use the current version at the time of writing which is Geant4.10.7.p2. We make use of the fact that the om-cmake function of om.bash is sensitive
to CMAKE_PREFIX_PATH envvar so that we can point to the directories where the libraries are installed and void having to rebuild them.  In principle just cut and paste the following line to a file change the envars of the different directories to match your system and source the resulting script.

    cd to the directory where you want to install Opticks (the WORK_DIR environmental variable will point to this directory). 
    
    
Here we are using a tagged snapshot of Opticks which can be found in github:

    git clone https://github.com/simoncblyth/opticks.git
    cd opticks
    git checkout tags/v0.1.6 -b v0.1.6-branch
    git status
    
The development version (a. k. a. the latest and greatest) can be found in the following repository 

    git clone https://bitbucket.org/simoncblyth/opticks.git
    
    

change opticks/optickscore/OpticksSwitches.h

so that:

    #define WITH_SKIPAHEAD 1

is set. 

    cat > setup_opticks.sh << +EOF
    # ----------------------------------------------------------------------------------------------------------------------
    # --- you need to modify the following environmental variables so that point to the specific directories on your system
    # --- 
    export WORK_DIR=/data2/wenzel/gputest_10.7.p02
    export OptiX_INSTALL_DIR=/home/wenzel/NVIDIA-OptiX-SDK-6.5.0-linux64
    export OPTICKS_COMPUTE_CAPABILITY=75
    export CUDA_INSTALL_DIR=/usr/local/cuda-11.3
    export CUDA_SAMPLES=${CUDA_INSTALL_DIR}/samples
    export G4INSTALL=/data2/wenzel/Geant4.10.07.p02_install
    . ${G4INSTALL}/bin/Geant4.sh
    export ROOTSYS=/data2/wenzel/root_install
    . ${ROOTSYS}/bin/thisroot.sh
    # ----------------------------------------------------------------------------------------------------------------------
    export LOCAL_BASE=${WORK_DIR}/local
    export CMAKE_PREFIX_PATH=${G4INSTALL}:${LOCAL_BASE}/opticks/externals:${OptiX_INSTALL_DIR}:${WORK_DIR}/opticks/cmake/Modules/:${WORK_DIR}/local/opticks:${WORK_DIR}/local/opticks:${WORK_DIR}/local/opticks/externals/
    export PYTHONPATH=$WORK_DIR
    export OPTICKS_HOME=${WORK_DIR}/opticks
    export PATH=${LOCAL_BASE}/bin:${PATH}
    export OPTICKS_PREFIX=${WORK_DIR}/local/opticks                            
    export OPTICKS_INSTALL_PREFIX=$LOCAL_BASE/opticks
    export OPTICKS_OPTIX_PREFIX=${OptiX_INSTALL_DIR}
    export OPTICKS_CUDA_PREFIX=${CUDA_INSTALL_DIR}
    export OPTICKS_EMBEDDED_COMMANDLINE_EXTRA="--rngmax 10 --rtx 1 --skipaheadstep 10000"
    opticks-(){ . ${OPTICKS_HOME}/opticks.bash && opticks-env $* ; }
    op(){ op.sh $* ; }
    o(){ cd $(opticks-home) ; hg st ; }
    # make sure to add the compiler options
    new=" -fPIC" 
    case ":${CXXFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) CXXFLAGS="$CXXFLAGS:$new"  ;;
    esac
    new=" -fPIC" 
    case ":${CFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) CFLAGS="$CFLAGS:$new"  ;;
    esac
    # speed up the make process
    new=" -j$(($(grep -c ^processor /proc/cpuinfo) - 1))" 
    case ":${MAKEFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) MAKEFLAGS="$MAKEFLAGS:$new"  ;;
    esac
    # deal with the $LD_LIBRARYPATH
    new=${OptiX_INSTALL_DIR}/lib64/
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${OPTICKS_HOME}/externals/lib
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${CUDA_INSTALL_DIR}/lib64/
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${LOCAL_BASE}/opticks/lib/
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    opticks-
    new=${CUDA_INSTALL_DIR}/bin
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${OPTICKS_HOME}/bin/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${OPTICKS_HOME}/ana/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${LOCAL_BASE}/opticks/lib/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${CUDA_SAMPLES}/bin/x86_64/linux/release/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    oinfo-(){
        echo 'LD_LIBRARY_PATH:';
        echo '================';
        echo  ${LD_LIBRARY_PATH}| tr : \\n;
        echo;
        echo 'PATH:';
        echo '=====';
        echo  ${PATH}| tr : \\n;
        echo;
        echo 'CMAKE_PREFIX_PATH:';
        echo '==================';
        echo  ${CMAKE_PREFIX_PATH}| tr : \\n;
        }
    dinfo-(){    
        nvidia-smi;
        ${CUDA_SAMPLES}/bin/x86_64/linux/release/deviceQuery
    }
    +EOF

    source setup_opticks.sh
    oinfo-
    dinfo-
    mkdir -p ${WORK_DIR}/local/opticks/externals/
    cd ${WORK_DIR}/local/opticks/externals/
    ln -s ${OptiX_INSTALL_DIR} OptiX
    cd ${WORK_DIR}
    opticks-externals-install >& install_ext.log &

scan the log file or any errors and correct them.

    cd ${WORK_DIR}
    opticks-full  >& install_full.log &

scan the log file or any errors and correct them. Before you run     opticks-t you want to to create the geocache using e.g. one of the gdml files provided by CaTS
(https://github.com/hanswenzel/CaTS)

    geocache-
    geocache-create- --gdmlpath  ${WORK_DIR}/local/opticks/opticksdata/export/juno1808/g4_00.gdml
    export OPTICKS_KEY=`output from above command"
    opticks-t
    
if the geocache-create- command doesn't work:


    ${WORK_DIR}/local/opticks/lib/OKX4Test --okx4test --g4codegen --deletegeocache --gdmlpath  ${WORK_DIR}/CaTS/gdml/simpleLArTPC.gdml
    export OPTICKS_KEY=`output from above command"
    opticks-t



# Building CaTS

Instructions for building and running CaTS can be found here:

https://github.com/hanswenzel/CaTS/blob/master/README.md
