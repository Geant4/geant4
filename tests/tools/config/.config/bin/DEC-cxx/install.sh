######################################
#
# g4system.U
#
G4SYSTEM="DEC-cxx"

#
# g4dirs.U
#
G4INSTALL="/usr/users/sadilov/geant4"
G4INCLUDE="/usr/users/sadilov/geant4/include"
G4BASE="/usr/users/sadilov/geant4/source"
G4WORKDIR="/usr/users/sadilov/geant4"
G4TMP="/usr/users/sadilov/geant4/tmp"
G4LIB="/usr/users/sadilov/geant4/lib"
G4BIN="/usr/users/sadilov/geant4/bin"
G4DATA="/usr/users/sadilov/geant4/data"
G4LEVELGAMMADATA="/usr/users/sadilov/geant4/data/PhotonEvaporation"

#
# g4clhep.U
#
CLHEP_BASE_DIR="/usr/clhep"
CLHEP_INCLUDE_DIR="/usr/local/include"
CLHEP_LIB_DIR="/usr/local/lib"
CLHEP_LIB="CLHEP"

#
# g4ospace
#
if [ Xn = Xy ] ; then
G4USE_OSPACE=1
else
G4USE_OSPACE=""
fi

G4OSPACE_BASE_DIR="/usr/users/sadilov/ObjectSpace"

#
# g4debug
#
G4DEBUG=""

#
# g4analysis
#
if [ Xn = Xy ] ; then
G4ANALYSIS_BUILD=1
else
G4ANALYSIS_BUILD=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_BUILD_LAB=1
else
G4ANALYSIS_BUILD_LAB=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_BUILD_JAS=1
else
G4ANALYSIS_BUILD_JAS=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_BUILD_LIZARD=1
else
G4ANALYSIS_BUILD_LIZARD=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_USE=1
else
G4ANALYSIS_USE=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_USE_LAB=1
else
G4ANALYSIS_USE_LAB=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_USE_JAS=1
else
G4ANALYSIS_USE_JAS=""
fi 

if [ Xn = Xy ] ; then
G4ANALYSIS_USE_LIZARD=1
else
G4ANALYSIS_USE_LIZARD=""
fi 

#
# g4ui
#
if [ Xy = Xy ] ; then
G4UI_BUILD_TERMINAL_SESSION=1
else
G4UI_BUILD_TERMINAL_SESSION=""
fi 

if [ Xy = Xy ] ; then
G4UI_USE_TERMINAL=1
else
G4UI_USE_TERMINAL=""
fi 

if [ Xn = Xy ] ; then
G4UI_BUILD_GAG_SESSION=1
else
G4UI_BUILD_GAG_SESSION=""
fi 

if [ Xn = Xy ] ; then
G4UI_USE_GAG=1
else
G4UI_USE_GAG=""
fi 

if [ Xn = Xy ] ; then
G4UI_BUILD_XAW_SESSION=1
else
G4UI_BUILD_XAW_SESSION=""
fi 

if [ Xn = Xy ] ; then
G4UI_USE_XAW=1
else
G4UI_USE_XAW=""
fi 

if [ Xn = Xy ] ; then
G4UI_BUILD_XM_SESSION=1
else
G4UI_BUILD_XM_SESSION=""
fi 

if [ Xn = Xy ] ; then
G4UI_USE_XM=1
else
G4UI_USE_XM=""
fi 

if [ Xn = Xy ] ; then
G4UI_BUILD_WIN32_SESSION=1
else
G4UI_BUILD_WIN32_session=""
fi 

if [ Xn = Xy ] ; then
G4UI_USE_WIN32=1
else
G4UI_USE_WIN32=""
fi 

if [ Xn = Xy ] ; then
G4UI_USE_TCSH=1
else
G4UI_USE_TCSH=""
fi 
#
# g4shared
#
G4LIB_BUILD_SHARED="n"
G4LIB_BUILD_STATIC="n"
G4LIB_USE_SHARED="n"
G4LIB_USE_STATIC="n"

if [ Xn = Xy ] ; then
G4LIB_BUILD_SHARED=1
else
G4LIB_BUILD_SHARED=""
fi 

if [ Xn = Xy ] ; then
G4LIB_BUILD_STATIC=1
else
G4LIB_BUILD_STATIC=""
fi 

if [ Xn = Xy ] ; then
G4LIB_USE_SHARED=1
else
G4LIB_USE_SHARED=""
fi 

if [ Xn = Xy ] ; then
G4LIB_USE_STATIC=1
else
G4LIB_USE_STATIC=""
fi 

#
# g4granular
#
if [ Xn = Xy ] ; then
G4LIB_USE_GRANULAR=1
else
G4LIB_USE_GRANULAR=""
fi 

#####################################################################
export G4SYSTEM

export G4INSTALL
export G4INCLUDE
export G4BASE
export G4WORKDIR
export G4TMP
export G4LIB
export G4BIN
export G4DATA
export G4LEVELGAMMADATA

export CLHEP_BASE_DIR
export CLHEP_INCLUDE_DIR
export CLHEP_LIB_DIR
export CLHEP_LIB

export G4USE_OSPACE
export G4OSPACE_BASE_DIR

export G4LIB_USE_GRANULAR

export G4DEBUG

export G4ANALYSIS_BUILD
export G4ANALYSIS_BUILD_LAB
export G4ANALYSIS_BUILD_JAS
export G4ANALYSIS_BUILD_LIZARD
export G4ANALYSIS_USE
export G4ANALYSIS_USE_LAB
export G4ANALYSIS_USE_JAS
export G4ANALYSIS_USE_LIZARD

export G4UI_BUILD_TERMINAL_SESSION
export G4UI_USE_TERMINAL
export G4UI_BUILD_GAG_SESSION
export G4UI_USE_GAG
export G4UI_BUILD_XAW_SESSION
export G4UI_USE_XAW
export G4UI_BUILD_XM_SESSION
export G4UI_USE_XM
export G4UI_BUILD_WIN32_SESSION
export G4UI_USE_WIN32
export G4UI_USE_TCSH

export G4LIB_BUILD_SHARED
export G4LIB_BUILD_STATIC
export G4LIB_USE_SHARED
export G4LIB_USE_STATIC

###############################################


#####################################################################
echo "On this machine the G4SYSTEM=$G4SYSTEM"

echo "On this machine the G4INSTALL=$G4INSTALL"
echo "On this machine the G4INCLUDE=$G4INCLUDE"
echo "On this machine the G4BASE=$G4BASE"
echo "On this machine the G4WORKDIR=$G4WORKDIR"
echo "On this machine the G4TMP=$G4TMP"
echo "On this machine the G4LIB=$G4LIB"
echo "On this machine the G4BIN=$G4BIN"
echo "On this machine the G4DATA=$G4DATA"
echo "On this machine the G4LEVELGAMMADATA=$G4LEVELGAMMADATA"

echo "On this machine the CLHEP_BASE_DIR=$G4CLHEP_BASE_DIR"
echo "On this machine the CLHEP_INCLUDE_DIR=$G4CLHEP_INCLUDE_DIR"
echo "On this machine the CLHEP_LIB_DIR=$G4CLHEP_LIB_DIR"
echo "On this machine the CLHEP_LIB=$G4CLHEP_LIB"

echo "On this machine the G4USE_OSPACE=$G4USE_OSPACE"
echo "On this machine the G4OSPACE_BASE_DIR=$G4OSPACE_BASE_DIR"

echo "On this machine the G4LIB_USE_GRANULAR=$G4LIB_USE_GRANULAR"

echo "On this machine the G4DEBUG=$G4DEBUG"

echo "On this machine the G4UI_BUILD_TERMINAL_SESSION=$G4UI_BUILD_TERMINAL_SESSION"
echo "On this machine the G4UI_USE_TERMINAL=$G4UI_USE_TERMINAL"
echo "On this machine the G4UI_BUILD_GAG_SESSION=$G4UI_BUILD_GAG_SESSION"
echo "On this machine the G4UI_USE_GAG=$G4UI_USE_GAG"
echo "On this machine the G4UI_BUILD_XAW_SESSION=$G4UI_BUILD_XAW_SESSION"
echo "On this machine the G4UI_USE_XAW=$G4UI_USE_XAW"
echo "On this machine the G4UI_BUILD_XM_SESSION=$G4UI_BUILD_XM_SESSION"
echo "On this machine the G4UI_USE_XM=$G4UI_USE_XM"
echo "On this machine the G4UI_BUILD_WIN32_SESSION=$G4UI_BUILD_WIN32_SESSION"
echo "On this machine the G4UI_USE_WIN32=$G4UI_USE_WIN32"
echo "On this machine the G4UI_USE_TCSH=$G4UI_USE_TCSH"

echo "On this machine the G4ANALYSIS_BUILD=$G4ANALYSIS_BUILD"
echo "On this machine the G4ANALYSIS_BUILD_LAB=$G4ANALYSIS_BUILD_LAB"
echo "On this machine the G4ANALYSIS_BUILD_JAS=$G4ANALYSIS_BUILD_JAS"
echo "On this machine the G4ANALYSIS_BUILD_LIZARD=$G4ANALYSIS_BUILD_LIZARD"
echo "On this machine the G4ANALYSIS_USE=$G4ANALYSIS_USE"
echo "On this machine the G4ANALYSIS_USE_LAB=$G4ANALYSIS_USE_LAB"
echo "On this machine the G4ANALYSIS_USE_JAS=$G4ANALYSIS_USE_JAS"
echo "On this machine the G4ANALYSIS_USE_LIZARD=$G4ANALYSIS_USE_LIZARD"

echo "On this machine the G4LIB_BUILD_SHARED=$G4LIB_BUILD_SHARED"
echo "On this machine the G4LIB_BUILD_STATIC=$G4LIB_BUILD_STATIC"
echo "On this machine the G4LIB_USE_SHARED=$G4LIB_USE_SHARED"
echo "On this machine the G4LIB_USE_STATIC=$G4LIB_USE_STATIC"

###############################################################
echo ""
echo "Starting installation..."
echo ""
cd /usr/users/sadilov/geant4/source
echo /usr/users/sadilov/work/g4conf/config/scripts
if [ Xn = Xy ] ; then
gmake global
if [ Xn = Xy ] ; then
gmake
fi
else
gmake
fi
if [ Xn = Xy ] ; then
gmake includes
fi
