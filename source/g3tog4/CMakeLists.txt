# - G3toG4 category build
include(sources.cmake)
geant4_add_category(G3toG4 MODULES G3toG4)

# OBSOLETE: build the rztog4 application
# It's not a Mac Bundle or WIN32 app (i.e. no GUI!!)
# Will need cernlib libraries, but no FindPackage for that yet...
# Hmm, means it'll probably be UNIX only...
# NB: Doesn't build on Linux due to undefined references - so do
# need CERNLIB
#
#include_directories(include)
#add_executable(rztog4
#    src/g3routines.F
#    src/jshape.F
#    src/rztog4.F
#    src/g3tog4.F
#    src/mztog4.F
#    src/tog4.F
#)
