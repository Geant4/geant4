// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#include "G4PreCompoundParameters.hh"


const G4double theLevelDensity = 0.125/MeV;

G4PreCompoundParameters G4PreCompoundParameters::thePreCompoundParameters;

G4PreCompoundParameters * G4PreCompoundParameters::GetAddress()
{ return &thePreCompoundParameters; }

