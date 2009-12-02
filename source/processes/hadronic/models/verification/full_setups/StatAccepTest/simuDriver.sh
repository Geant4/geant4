#!/bin/sh
#
#----------------------------------------------------------------------------
# Last update: 02-Dec-2009.
#
# This Bash shell script has the following 11 parameters:
#
#   1)  Geant4 reference; e.g. 4.6.2.ref03 , or a local directory. 
#   2)  Flag to know whether the first reference should be simulated or not.
#   3)  Geant4 reference; e.g. 4.6.2.ref04 , or empty, or a local directory.
#   4)  Flag to know whether the second reference should be simulated or not.
#   5)  Flag to know whether the statistical tests should be run or not.
#   6)  Physics List; e.g. LHEP
#   7)  Calorimeter type; e.g. FeSci
#   8)  Particle type; e.g. pi+
#   9)  Beam Energy; e.g. 20GeV
#   10) Number of Events; e.g. 5k
#   10) Bfield; e.g. 4.0tesla
#
# This script invokes the Python  build.py  that writes the setup
# and the Geant4 command file, and then run the simulation(s).
# After that, this script eventually (if the flag is on) calls the 
# Python script  driveStatTest.py  that does the Statistical tests.
#
# Notice that the "LABEL" which specifies all the input parameters
# does not contain the Bfield in the default case of zero magnetic field.
#
#----------------------------------------------------------------------------
#
echo ' ========== START simuDriver.sh ========== '
#
export REF1=$1
export SIM_REF1=$2
export REF2=$3
export SIM_REF2=$4
export RUN_STAT=$5
export PHYSICS=$6
export CALORIMETER=$7
export PARTICLE=$8
export ENERGY=$9
export EVENTS=${10}
export BFIELD=${11}
#
echo ' REF1        =' $REF1
echo ' SIM_REF1    =' $SIM_REF1
echo ' REF2        =' $REF2
echo ' SIM_REF2    =' $SIM_REF2
echo ' RUN_STAT    =' $RUN_STAT
echo ' PHYSICS     =' $PHYSICS
echo ' CALORIMETER =' $CALORIMETER
echo ' PARTICLE    =' $PARTICLE
echo ' ENERGY      =' $ENERGY
echo ' EVENTS      =' $EVENTS
echo ' BFIELD      =' $BFIELD
#
export UNAMEI=`uname -i`
echo ' UNAMEI      =' $UNAMEI

#--- Run the first reference ---
#
( if [ X$SIM_REF1 == XYes ] ; then
#
    ###echo ' I AM HERE 1 ' ;
    export REF=$REF1 ;
    export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
    if [ X$BFIELD != X ] ; then
      if [ X$BFIELD != X0 ] ; then
        if [ X$BFIELD != X0. ] ; then
          if [ X$BFIELD != X0.0 ]; then
	      export LABEL=$LABEL-B$BFIELD ;
          fi
        fi
      fi
    fi
    ###echo ' 1) LABEL=' $LABEL
    python build.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS $BFIELD $UNAMEI ;
    if [ $? != 0 ] ; then
	echo ' ***ERROR*** from: python build.py ... ! exitCode = 11' ;
	exit 11 ;
    fi
    mv run.g4 run.g4-$LABEL ;
    mv setup.sh setup.sh-$LABEL ;
    . setup.sh-$LABEL ;
    if [ $? != 0 ] ; then
	echo ' ***ERROR*** from: . setup.sh-... ! exitCode = 12' ; 
	exit 12 ;	
    fi
#
    echo '  ' ;
    echo '--- Check platform / environment --- ' ;
    echo '*** g++ -v ***' ;                    g++ -v                    ; echo ' ' ;
    echo '*** which g++ ***' ;                 which g++                 ; echo ' ' ;
    echo '*** uname -a ***' ;                  uname -a                  ; echo ' ' ;
    echo '*** cat /etc/issue ***' ;            cat /etc/issue ;        
    echo '*** cat /etc/cpuinfo ***' ;          cat /proc/cpuinfo ;
    echo '*** DIR_INSTALLATIONS = '            $DIR_INSTALLATIONS        ; echo ' ' ;
    echo '*** ls -lh $DIR_INSTALLATIONS ***' ; ls -lh $DIR_INSTALLATIONS ; echo ' ' ;
    echo '*** G4INSTALL = '                    $G4INSTALL                ; echo ' ' ;
    echo '*** ls -lh $G4INSTALL ***' ;         ls -lh $G4INSTALL         ; echo ' ' ;
    echo '*** PWD = '                          $PWD                      ; echo ' ' ;
    echo '*** ls -lh $PWD *** ' ;              ls -lh $PWD ;
    echo '------------------------------------ ' ;    
    echo '  ' ;
#
    echo '  '; echo ' G4INSTALL = ' $G4INSTALL; echo ' running REF = ' $REF ; echo '  ' ;

    mainStatAccepTest-$REF-$RELEASENAME-`uname -i` run.g4-$LABEL > output.log-$LABEL 2>&1 ;

    EXITCODE=$? ;
    if [ $EXITCODE != 0 ] ; then
	echo ' ***ERROR*** from: mainStatAccepTest-... run.g4-... ! exitCode = 14' ;  
	rm -rf tmp/ ;
	exit $EXITCODE ;
    fi
    mv ntuple.root ntuple.root-$LABEL ;
#
    echo ' ' ;
    echo '--- Check results after running 1st reference ---' ; 
    echo '*** ls -lth ***' ;  ls -lth ;
    echo '-------------------------------------------------' ;
    echo ' ' ;
#
  fi )
#
if [ $? != 0 ] ; then
    isFirstBad=Yes ;
    ###echo ' isFirstBad = ' $isFirstBad ;
fi
#
#--- Run the second reference ---
#
( if [ X$SIM_REF2 == XYes ] ; then
#
    ###echo ' I AM HERE 2 ' ;
    export REF=$REF2 ; 
    export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
    if [ X$BFIELD != X ] ; then
      if [ X$BFIELD != X0 ] ; then
        if [ X$BFIELD != X0. ] ; then
          if [ X$BFIELD != X0.0 ]; then
	      export LABEL=$LABEL-B$BFIELD ;
          fi
        fi
      fi
    fi
    ###echo ' 2) LABEL=' $LABEL ;
    python build.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS $BFIELD $UNAMEI ;
    if [ $? != 0 ] ; then
	echo ' ***ERROR*** from: python build.py ... ! exitCode = 21' ;
	exit 21	;
    fi
    mv run.g4 run.g4-$LABEL ;
    mv setup.sh setup.sh-$LABEL ;
    . setup.sh-$LABEL ;
    if [ $? != 0 ] ; then
	echo ' ***ERROR*** from: . setup.sh-... ! exitCode = 22' ; 
	exit 22 ;	
    fi
    echo '  '; echo ' G4INSTALL = ' $G4INSTALL; echo ' running REF = ' $REF ; echo '  ' ;

    mainStatAccepTest-$REF-$RELEASENAME-`uname -i` run.g4-$LABEL > output.log-$LABEL 2>&1 ;

    EXITCODE=$? ;
    if [ $EXITCODE != 0 ] ; then
	echo ' ***ERROR*** from: mainStatAccepTest-... run.g4-... ! exitCode = 24' ;   
	rm -rf tmp/ ;
	exit $EXITCODE ;
    fi
    mv ntuple.root ntuple.root-$LABEL ;
#
    echo ' ' ;
    echo '--- Check results after running 2nd reference ---' ; 
    echo '*** ls -lth ***' ;  ls -lth ;
    echo '-------------------------------------------------' ;
    echo ' ' ;
#
  fi )
#
if [ $? != 0 ] ; then
    isSecondBad=Yes ;
    ###echo ' isSecondBad = ' $isSecondBad ;
fi
#
#--- Run the statistical tests ---
#
( if [ X$RUN_STAT == XYes ] ; then
#
    if [[ X$isFirstBad != XYes && X$isSecondBad != XYes ]] ; then 

        ###echo ' I AM HERE 3 ' ;
	export LABEL=$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
	if [ X$BFIELD != X ] ; then
	    if [ X$BFIELD != X0 ] ; then
		if [ X$BFIELD != X0. ] ; then
		    if [ X$BFIELD != X0.0 ]; then
			export LABEL=$LABEL-B$BFIELD ;
		    fi
		fi
	    fi
	fi
        ###echo ' 3) LABEL=' $LABEL ;
	. setup.sh-$REF1-$LABEL ;
	if [ $? != 0 ] ; then
	    echo ' ***ERROR*** from: . setup.sh-... ! exitCode = 31' ;
	    exit 31 ;	
	fi

	python driveStatTest.py $REF1 $REF2 $LABEL ; 
	if [ $? != 0 ] ; then
	    echo ' ***ERROR*** from: python driveStatTest.py ... ! exitCode = 33' ;   
	    exit 33 ;
	fi
#
	echo ' ' ;
	echo '--- Check results after running statistical test ---' ;
	echo '*** ls -lth ***' ;  ls -lth ;
	echo '-------------------------------------------------' ;
	echo ' ' ;
#
    fi

  fi )
#
if [ $? != 0 ] ; then
    isThirdBad=Yes ;
    ###echo ' isThirdBad = ' $isThirdBad ;
fi
#
echo ' ========== END simuDriver.sh ========== '
#
if [[ X$isFirstBad == XYes || X$isSecondBad == XYes || X$isThirdBad == XYes ]] ; then 
    echo ' ************************************************** ' ;
    echo ' ******       simuDriver.sh  FAILED!         ****** ' ;
    echo ' isFirstBad  = ' $isFirstBad ;
    echo ' isSecondBad = ' $isSecondBad ;
    echo ' isThirdBad  = ' $isThirdBad ;
    echo ' ************************************************** ' ;
    exit 69 ;
else
    echo ' ************************************************** ' ;
    echo ' ******       simuDriver.sh  OK!             ****** ' ;
    echo ' ************************************************** ' ;
fi
