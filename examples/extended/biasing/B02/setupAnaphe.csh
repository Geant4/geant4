#!/usr/bin/tcsh -x

#set verbose = 0

set release = "5.0.1"
set redHat72 = "redhat72"
set redHat61 = "redhat61"
set redhat73 = "redhat73"
set gcc295 = "gcc-2.95.2"
set gcc296 = "gcc-2.96"
set gcc32  = "gcc-3.2"
set platform = "$redHat72/$gcc295" # default
set COMP = `g++ --version | tail -1`
set os = `uname`
set osVersion = `uname -r`
set unavailable = "Sorry, Anaphe $release is not available for $os, compiler $COMP"



# deal with non-Linux (= Solaris CC 5.3) first
if ($os != "Linux" && ) then
  if ($os == "SunOS" && $osVersion == "5.8") then
    set platform = "sun4x_58/CC-5.3"
  else
    echo $unavailable
    exit
  endif
endif


# Linux next. Find Redhat/SuSE flavours
set is61 = "no"
if { grep "6.1"  /etc/issue > /dev/null } then
  set is61 = "yes"
endif
set is72 = "no"
if { grep "7.2"  /etc/issue > /dev/null } then
  set is72 = "yes"
endif
set is73 = "no"
if { grep "7.3"  /etc/issue > /dev/null } then
  set is73 = "yes"
endif
set isSuSE = "no"
if { grep "SuSE" /etc/issue > /dev/null } then
  set isSuSE = "yes"
endif
# Constraint: all these vars are now "yes" or "no"


# deal with SuSE
if ($isSuSE == "yes") then
  echo "Found SuSE distribution"
  if ($is72 == "yes") then  # SuSE 7.2 is compatible with RH61
    # use Red Hat 6.1 / gcc-2.95.2  for SuSE 7.2 (gcc-2.95.3)
    set platform = "$redHat61/$gcc295"
  else
    echo "Unknown version of SuSE - please report to HepLib.Support@cern.ch"
    echo "/etc/issue contains:"
    cat /etc/issue
    exit 
  endif
endif


# deal with non-SuSE - assumed to be Red Hat (but might be wrong..!)
if ($isSuSE == "no") then
  # first deal with Red Hat 6.1
  if ($is61 == "yes") then
    if ($COMP == "gcc-2.96") then
      echo $unavailable
      exit
    else if ($COMP == "2.95.2") then
      set platform = "redhat61/gcc-2.95.2"
    else
      set platform = "redhat61/egcs_1.1.2" # assumption!
    endif
  endif
  # now deal with Red Hat 7.2
  if ($is72 == "yes") then
    if ($COMP == "gcc-2.96") then
      echo $unavailable
      exit
    else if ($COMP == "2.95.2") then
      set platform = "redhat72/gcc-2.95.2"
    else
      set platform = "redhat72/gcc-3.2" # assumption!
    endif
  endif
  # finally deal with Red Hat 7.3
  else if ($is73 == "yes") then
    if ( $COMP == "gcc-2.96" ) then
      echo $unavailable
      exit
    else if ( $COMP == "gcc-2.95.2" ) then
      set platform = "redhat73/gcc-2.95.2"
    else
      set platform = "redhat73/gcc-3.2" # assumption!
    endif
  endif
else
  echo $unavailable
endif



setenv PLATF $platform



if ($?AIDA_DIR == 0) then
   setenv AIDA_DIR /afs/cern.ch/sw/contrib/AIDA/3.0/src/cpp
endif

if ($?ANAPHETOP == 0) then
   setenv ANAPHETOP /afs/cern.ch/sw/lhcxx
endif

# add to PATH and LD_LIBRARY_PATH, but only if not already present ...
set pubDom = "${ANAPHETOP}/specific/${PLATF}/PublicDomainPackages/2.0.0"

if ( { echo $PATH | grep $pubDom } == 0 ) then
  setenv PATH "${pubDom}/bin:${PATH}"
endif

if ( $?LD_LIBRARY_PATH == 0 ) then
  setenv LD_LIBRARY_PATH "${pubDom}/lib"
else if ( { echo $LD_LIBRARY_PATH | grep $pubDom } == 0 ) then
  setenv LD_LIBRARY_PATH "${pubDom}/lib:${LD_LIBRARY_PATH}"
endif

set releaseDirectory = "${ANAPHETOP}/specific/${PLATF}/${release}"

if ( { echo $PATH | grep $releaseDirectory } == 0 ) then
  setenv PATH "${releaseDirectory}/bin:${PATH}"
endif

if ( { echo $LD_LIBRARY_PATH | grep $releaseDirectory } == 0 ) then
  setenv LD_LIBRARY_PATH "${releaseDirectory}/lib:${LD_LIBRARY_PATH}"
endif

if ( $?GRACE_HOME == 0 ) then
  setenv GRACE_HOME "${pubDom}/grace"
endif

if ( { echo $PATH | grep $GRACE_HOME } == 0 ) then
  setenv PATH "${GRACE_HOME}/bin:${PATH}"
endif
