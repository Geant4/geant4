//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StableIsotopes.hh,v 1.6 2001-08-01 17:12:42 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4StableIsotopes_h
#define G4StableIsotopes_h 1

// Class Description
// class utility that knows about all naturally occuring isotopes.; 
// to be used in your process implementation in case you need this.
// Class Description - End


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

private:

static const G4int protonCount[92];
static const G4String elementName[92];
static const G4int nIsotopes[92];
static const G4int start[92];
static const G4int nucleonCount[287];
static const G4double abundance[287];
};

#endif
