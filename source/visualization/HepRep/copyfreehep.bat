set HEPREP=%FREEHEP%\hep\graphics\heprep
copy %HEPREP%\org\freehep\jni\FreeHepTypes.h include
copy %HEPREP%\geant4\include\*.hh include
copy %HEPREP%\streamer-c++\*.h include
copy %HEPREP%\include\HEPREP\*.h include\HEPREP
copy %HEPREP%\geant4\*.cc src
copy %HEPREP%\streamer-c++\*.cc src

