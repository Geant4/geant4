set HEPREP=%FREEHEP%\hep\graphics\heprep
copy "%FREEHEP%\org\freehep\jni\FreeHepTypes.h"     "include\ref"
REM copy "%FREEHEP%\org\freehep\zip\zipios++\*.h"       "include\zipios++"
REM copy "%FREEHEP%\org\freehep\zip\*.h"                "include\zipios++"
REM copy "%FREEHEP%\org\freehep\zip\*.cpp"              "src\*.cc"
copy "%HEPREP%\c++\include\*.h"                     "include\ref"
copy "%HEPREP%\include\HEPREP\*.h"                  "include\HEPREP"
copy "%HEPREP%\c++\*.cc"                            "src"
del src\MultiWriteTest.cc
del src\HepRepTest.cc
del src\HepRepExample.cc
REM del src\test_zipoutputstream.cc
REM del src\test_gzipoutputstream.cc
