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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include <iomanip.h>

#include "globals.hh"
#include "Randomize.hh"

#include "G4CascadSpecialFunctions.hh"

#include "vector"

G4int crossSections();

int main(int argc, char **argv ) {
  crossSections();   

  return 0;       
};

G4int crossSections() {   // print cross section data


  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> crossSections() " << G4endl;
  }

  if (verboseLevel > 1) {
    G4cout << " MeV: " << MeV << " GeV: " << GeV << G4endl;
  }

  // 100  <> 10 GeV
  const G4int types[] = { 1, 2, 3, 5, 7};

  for (G4int iE = 1; iE < 151; iE++) { 

    if (verboseLevel > 0) {
      cout.precision(4);
      G4double e = G4double(iE) / 10.0;
      G4cout << setw(9)  << e;

	for (G4int j = 0; j < 5; j++) {
	  G4cout << setw(9)  << G4CascadSpecialFunctions::crossSection(e, types[j]);
	}}

    G4cout << G4endl;
  }

  return 0;
};






