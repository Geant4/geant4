#
# Setup for Unix (Linux, Solaris, MacOS X)
#

#
# Change this to point to your Java installation
#
#export JDK_HOME="/usr/java/jdk1.4.1"
#
#
# No User Configuration below.
#
if test ! "${JDK_HOME:+x}"
then
    export JDK_HOME=/usr/java/jdk1.4.1
fi
if ! test -e $JDK_HOME/bin/java
then
   echo "You need to install Java and/or modify script to set JDK_HOME correctly"
   return
fi

if test ! "${GEANT4_HOME:+x}"
then
if test ! "${G4INSTALL:+x}"
then
    GEANT4_HOME="../../../.."
else
    GEANT4_HOME=$G4INSTALL
fi
fi
export JLIBPATH=$GEANT4_HOME/extlib/java
if ! test -e $JLIBPATH/aida.jar
then
   echo "The AIDA.jar files seem to be missing, make sure JAIDA-Geant4 patch has been applied"
	echo "and that GEANT4_HOME is set correctly"
   return
fi

export G4ANALYSIS_USE=1
export G4ANALYSIS_AIDA_CONFIG_LIBS="-L$GEANT4_HOME/extlib/$G4SYSTEM -lAIDAJNI -lFHJNI"
touch src/DummyAnalysisFactory.cc

export JVM_ARGS="-Dhep.aida.IAnalysisFactory=jas.aida.gui.JASGUIAnalysisFactory"
case $G4SYSTEM in
Linux*)
    LD_LIBRARY_PATH=$JDK_HOME/jre/lib/i386:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$JDK_HOME/jre/lib/i386/client:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    G4ANALYSIS_AIDA_CONFIG_LIBS="$G4ANALYSIS_AIDA_CONFIG_LIBS -L$JDK_HOME/jre/lib/i386/client -ljvm"
    ;;
SUN*)
    LD_LIBRARY_PATH=$JDK_HOME/jre/lib/sparc:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$JDK_HOME/jre/lib/sparc/client:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    G4ANALYSIS_AIDA_CONFIG_LIBS="$G4ANALYSIS_AIDA_CONFIG_LIBS -L$JDK_HOME/jre/lib/sparc -ljvm"
    ;;
Darwin*)
    G4ANALYSIS_AIDA_CONFIG_LIBS="$G4ANALYSIS_AIDA_CONFIG_LIBS -framework JavaVM"
    ;;
*)
    break
esac

export CLASSPATH="$JLIBPATH/JASAIDA.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/aida.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/aida-dev.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/freehep-hep.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/freehep-base.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/openide-lookup.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/jas.jar"
CLASSPATH="$CLASSPATH:$JLIBPATH/jel.jar"

echo "JDK_HOME set to $JDK_HOME"
echo "JVM_ARGS set to $JVM_ARGS"
echo "CLASSPATH set to $CLASSPATH"
if test "${GEANT4_HOME:+x}"
then
    echo "LD_LIBRARY_PATH set to $LD_LIBRARY_PATH"
fi
echo "G4ANALYSIS_AIDA_CONFIG_LIBS set to $G4ANALYSIS_AIDA_CONFIG_LIBS"
