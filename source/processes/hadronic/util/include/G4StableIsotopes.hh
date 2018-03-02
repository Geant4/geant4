//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
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

G4String GetName(G4int Z);
G4int GetNumberOfIsotopes(G4int Z);
G4int GetFirstIsotope(G4int Z);
G4int GetIsotopeNucleonCount(G4int number);
G4double GetAbundance(G4int number);
G4int GetProtonCount(G4int Z);

private:

static const G4int protonCount[92];
static const G4String elementName[92];
static const G4int nIsotopes[92];
static const G4int start[92];
static const G4int nucleonCount[287];
static const G4double abundance[287];
};

#endif
