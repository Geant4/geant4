#!/bin/tcsh
cd $G4BUILD
make G4dna || cd - && exit
cd -
make $1
