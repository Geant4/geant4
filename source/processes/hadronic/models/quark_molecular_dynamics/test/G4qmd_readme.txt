#  Implementation of quark molecular dynamics (qMD) into GEANT4
#  ------------------------------------------------------------
#
#
#  2000-08-04, Stefan Scherer
#

G4qmdTest.cc reads in qMD data file, 
initializes qMD (by instantiating a class G4qmd)
and runs the code - 
output is writte to STDOUT, in qMD ascii output format.

WARNING: code produces runtime errors...


for possible input files, see *.dat files in

   /afs/cern.ch/user/s/sscherer/public/qmd/data

the meaning of the entries in the qMD output is explained in

   /afs/cern.ch/user/s/sscherer/public/qmd/qmd_format.txt

for further physics information, see links at
   
  http://cern.ch/Stefan.Scherer




