#
# Change this to point to your Java installation
#
export JDK_HOME=/usr/java/jdk1.4.1
if ! test -e $JDK_HOME/bin/java
then
   echo "You need to install Java and/or modify script to set JDK_HOME correctly"
   return
fi
export G4ANALYSIS_USE=1
export G4ANALYSIS_AIDA_CONFIG_LIBS="-lFREEHEP -L$JDK_HOME/jre/lib/i386/client -ljvm"
touch src/DummyAnalysisFactory.cc

export JVM_ARGS="-Dhep.aida.IAnalysisFactory=jas.aida.gui.JASGUIAnalysisFactory"
export LD_LIBRARY_PATH=$JDK_HOME/jre/lib/i386/client:$JDK_HOME/jre/lib/i386

export JLIBPATH=$G4INSTALL/lib/java
if ! test -e $JLIBPATH/aida.jar
then
   echo "The AIDA.jar files seem to be missing, make sure JAIDA-G4 patch has been applied"
	echo "and that G4INSTALL is set correctly"
   return
fi
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
echo "LD_LIBRARY_PATH set to $LD_LIBRARY_PATH"
echo "G4ANALYSIS_AIDA_CONFIG_LIBS set to $G4ANALYSIS_AIDA_CONFIG_LIBS"
