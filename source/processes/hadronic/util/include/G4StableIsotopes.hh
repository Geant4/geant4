// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StableIsotopes.hh,v 1.1 1999-01-07 16:13:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4StableIsotopes_h
#define G4StableIsotopes_h 1

// class, that knows the stable isotopes for all elements.
// H.P. Wellisch, 21. Nov. 1997
// Accessable by name, and Z
//
// To loop over all isotopes for element with protonount Z
// G4StableIsotopes theIso;
// for (G4int i=0; i<theIso.GetNumberOfIsotopes(); i++)
// {
//    G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(Z)+i);
// }
//
#include "globals.hh"

class G4StableIsotopes
{

public:

G4String GetName(G4int Z) {return elementName[Z-1];}
G4int GetNumberOfIsotopes(G4int Z) {return nIsotopes[Z-1];}
G4int GetFirstIsotope(G4int Z) {return start[Z-1];}
G4int GetIsotopeNucleonCount(G4int number) {return nucleonCount[number];}
G4double GetAbundance(G4int number) {return abundance[number];}
G4int GetProtonCount(G4int Z) {return protonCount[Z-1];}

public:

static const G4int protonCount[92];
static const G4String elementName[92];
static const G4int nIsotopes[92];
static const G4int start[92];
static const G4int nucleonCount[287];
static const G4double abundance[287];
};

#endif
