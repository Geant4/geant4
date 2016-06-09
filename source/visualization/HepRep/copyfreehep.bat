set HEPREP=%FREEHEP%\hep\graphics\heprep
copy "%FREEHEP%\org\freehep\jni\FreeHepTypes.h"     "include\ref"
copy "%FREEHEP%\org\freehep\zip\zipios++\*.h"       "include\zipios++"
copy "%FREEHEP%\org\freehep\zip\*.h"                "include\zipios++"
copy "%FREEHEP%\org\freehep\zip\*.cpp"              "src\*.cc"
copy "%HEPREP%\geant4\include\*.hh"                 "include"
copy "%HEPREP%\c++\include\*.h"                     "include\ref"
copy "%HEPREP%\include\HEPREP\*.h"                  "include\HEPREP"
copy "%HEPREP%\geant4\*.cc"                         "src"
copy "%HEPREP%\c++\*.cc"                            "src"
del src\MultiWriteTest.cc
del src\HepRepTest.cc
del src\HepRepExample.cc
del src\test_zipoutputstream.cc
del src\test_gzipoutputstream.cc
