#
# Setup for Unix (Linux, Solaris, MacOS X)
#

#
# Change this to point to your Java installation
#
#setenv JDK_HOME "/usr/java/jdk1.4.1"
#
#
# No User Configuration below.
#
if (! $?JDK_HOME) then
   setenv JDK_HOME /usr/java/jdk1.4.1
endif
if (! -e $JDK_HOME/bin/java) then
   echo "You need to install Java and/or modify script to set JDK_HOME correctly"
   goto end
endif

if (! $?GEANT4_HOME ) then
if (! $?G4INSTALL ) then
    setenv GEANT4_HOME "../../../.."
else
    setenv GEANT4_HOME $G4INSTALL
endif
endif
setenv JLIBPATH "$GEANT4_HOME/extlib/java"
if (! -e $JLIBPATH/aida.jar) then
   echo "The AIDA.jar files seem to be missing, make sure JAIDA-Geant4 patch has been applied"
	echo "and that GEANT4_HOME is set correctly"
   goto end
endif

setenv G4ANALYSIS_USE 1
setenv G4ANALYSIS_AIDA_CONFIG_LIBS "-L$GEANT4_HOME/extlib/$G4SYSTEM -lAIDAJNI -lFHJNI"
touch src/DummyAnalysisFactory.cc

setenv JVM_ARGS "-Dhep.aida.IAnalysisFactory=jas.aida.gui.JASGUIAnalysisFactory"
switch ($G4SYSTEM)
    case "Linux*":
        setenv LD_LIBRARY_PATH $JDK_HOME/jre/lib/i386:${LD_LIBRARY_PATH}
        setenv LD_LIBRARY_PATH $JDK_HOME/jre/lib/i386/client:${LD_LIBRARY_PATH}
        setenv G4ANALYSIS_AIDA_CONFIG_LIBS "$G4ANALYSIS_AIDA_CONFIG_LIBS -L$JDK_HOME/jre/lib/i386/client -ljvm"
        breaksw
    case "SUN*":
        setenv LD_LIBRARY_PATH $JDK_HOME/jre/lib/sparc:${LD_LIBRARY_PATH}
        setenv LD_LIBRARY_PATH $JDK_HOME/jre/lib/sparc/client:${LD_LIBRARY_PATH}
        setenv G4ANALYSIS_AIDA_CONFIG_LIBS "$G4ANALYSIS_AIDA_CONFIG_LIBS -L$JDK_HOME/jre/lib/sparc -ljvm"
        breaksw
    case "Darwin*":
        setenv G4ANALYSIS_AIDA_CONFIG_LIBS "$G4ANALYSIS_AIDA_CONFIG_LIBS -framework JavaVM"
        breaksw
    default:
        breaksw
endsw

setenv CLASSPATH "$JLIBPATH/JASAIDA.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/aida.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/aida-dev.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/freehep-hep.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/freehep-base.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/openide-lookup.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/jas.jar"
setenv CLASSPATH "${CLASSPATH}:$JLIBPATH/jel.jar"

echo "JDK_HOME set to $JDK_HOME"
echo "JVM_ARGS set to $JVM_ARGS"
echo "CLASSPATH set to $CLASSPATH"
if ( $?LD_LIBRARY_PATH ) then
    echo "LD_LIBRARY_PATH set to $LD_LIBRARY_PATH"
endif
echo "G4ANALYSIS_AIDA_CONFIG_LIBS set to $G4ANALYSIS_AIDA_CONFIG_LIBS"
end:
