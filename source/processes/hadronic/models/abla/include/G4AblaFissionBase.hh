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
// ABLAXX statistical de-excitation model
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4AblaFissionBase_hh
#define G4AblaFissionBase_hh 1

#ifdef ABLAXX_IN_GEANT4_MODE
#include "globals.hh"
#else
#include "G4INCLGeant4Compat.hh"
#endif
//#include "G4InclUtils.hh"

/*
 * Abstract interface to fission models.
 */

class G4AblaFissionBase {

public:
  G4AblaFissionBase();
  virtual ~G4AblaFissionBase();

  virtual void doFission(G4double &A, G4double &Z, G4double &E,
			 G4double &A1, G4double &Z1, G4double &E1, G4double &K1,
			 G4double &A2, G4double &Z2, G4double &E2, G4double &K2) = 0;

  void setVerboseLevel(G4int level) {
    verboseLevel = level;
  }

  void about() {
    //    G4cout << ";; " << aboutModel << G4endl;
  }

  void setAboutString(G4String anAbout) {
    aboutModel = anAbout;
  }

private:
  G4int verboseLevel;
  G4String aboutModel;
};

#endif
