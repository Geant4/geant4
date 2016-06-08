#
# Change this to point to your Java installation
#
setenv JDK_HOME /usr/java/jdk1.4.1
if (! -e $JDK_HOME/bin/java) then
   echo "You need to install Java and/or modify script to set JDK_HOME correctly"
   goto end
endif
setenv G4ANALYSIS_USE 1
setenv G4ANALYSIS_AIDA_CONFIG_LIBS "-lAIDAJNI -lFHJNI -L$JDK_HOME/jre/lib/i386/client -ljvm"
touch src/DummyAnalysisFactory.cc

setenv JVM_ARGS "-Dhep.aida.IAnalysisFactory=jas.aida.gui.JASGUIAnalysisFactory"
setenv LD_LIBRARY_PATH $JDK_HOME/jre/lib/i386/client:$JDK_HOME/jre/lib/i386

setenv JLIBPATH $G4INSTALL/lib/java
if (! -e $JLIBPATH/aida.jar) then
   echo "The AIDA.jar files seem to be missing, make sure JAIDA-G4 patch has been applied"
	echo "and that G4INSTALL is set correctly"
   goto end
endif
setenv CLASSPATH "$JLIBPATH/JASAIDA.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/aida.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/aida-dev.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/freehep-hep.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/freehep-base.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/openide-lookup.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/jas.jar"
setenv CLASSPATH "$CLASSPATH:$JLIBPATH/jel.jar"

echo "JDK_HOME set to $JDK_HOME"
echo "JVM_ARGS set to $JVM_ARGS"
echo "CLASSPATH set to $CLASSPATH"
echo "LD_LIBRARY_PATH set to $LD_LIBRARY_PATH"
echo "G4ANALYSIS_AIDA_CONFIG_LIBS set to $G4ANALYSIS_AIDA_CONFIG_LIBS"
end:
