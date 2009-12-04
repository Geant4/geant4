#!/bin/tcsh

set macroName = $1
set startingSeed = $2

set macroFile = macros/${macroName}.mac
set outputFile = output/${macroName}_$startingSeed

${G4WORKDIR}/bin/${G4SYSTEM}/electronScattering2 $macroFile $startingSeed $outputFile

exit
