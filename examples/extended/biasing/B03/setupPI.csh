# #Id: setupPI.csh,v 1.6 2004/06/19 17:22:46 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-08-00 $
# -------------------------------------------------------------------

# before sourcing this script make sure you have set the
# environment variables according to the description in README
# NOTE:: B03 cannot work with PI Ultralite, because it needs the Python include files which are not comprised
# in the Ultralite version

# setting in case PI and analysis is used
if (${?PI_BASE_DIR} == 1) then
   setenv PI_VERSION 1_3_0
   setenv PI_VER 1.3.0
   setenv PLATF slc3_ia32_gcc323
   setenv PYTHONVERSION 2.3
   setenv PI_DIR ${PI_BASE_DIR}/${PI_VER}/app/releases/PI/PI_${PI_VERSION}/${PLATF}
   setenv PATH ${PATH}:${PI_DIR}/bin
   setenv PATH ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PLATF}/bin:$PATH
   eval `aida-config -r csh`
#
# python from PI
#
  setenv PYTHON_INCLUDE_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PLATF}/include/python${PYTHONVERSION}
  setenv PYTHON_LIB_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PLATF}/lib/python${PYTHONVERSION}/config
#
# python for PI
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
  setenv PYTHONPATH ${PI_DIR}/python:$PYTHONPATH
  setenv PYTHONPATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_6_0/${PLATF}/python:$PYTHONPATH
#
# SWIG settings
setenv SWIGVERSION 1.3.15
#
else   # not using PI

  echo "-- WARNING: histograms are not activated !"
  echo "            Either PI_BASE_DIR are not set."
  #
  if (   ${?PYTHONVERSION} == 1 && \
         ${?PYTHON_BASE_DIR} == 1 && \
         ${?SWIG_BASE_DIR} == 1 && \
         ${?SWIG_VERSION} == 1 && \
         ${?CLHEP_BASE_DIR} == 1) then  # settings without PI
    setenv PYTHON_LIB_DIR ${PYTHON_BASE_DIR}/lib/python${PYTHONVERSION}/config
    setenv PYTHON_INCLUDE_DIR ${PYTHON_BASE_DIR}/include/python${PYTHONVERSION}
    #
    if ( ! -d $PYTHON_LIB_DIR ) then
      echo -- ERROR: no pyhton lib/config directory: $PYTHON_LIB_DIR
    endif
    if ( ! -d $PYTHON_INCLUDE_DIR ) then
      echo -- ERROR: no pyhton include directory: $PYTHON_INCLUDE_DIR
    endif
  else    # environment not completed in case no PI is used
    echo "-- ERROR: PYTHONVERSION or PYTHON_BASE_DIR or CLHEP_BASE_DIR or SWIG_VERSION or SWIG_BASE_DIR not set!"    
    exit
  endif

endif  # end settings for not using PI
#
if ( ${?LD_LIBRARY_PATH} == 0 ) then
  setenv LD_LIBRARY_PATH
endif
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/exampleB03
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PLATF}/lib:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PI_BASE_DIR}/${PI_VER}/external/Boost/1.31.0/${PLATF}/lib
setenv B03_BASE `pwd`
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
setenv PYTHONPATH ${PYTHONPATH}:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${PYTHONPATH}:${B03_BASE}
