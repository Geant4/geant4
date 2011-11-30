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
#ifndef G4InclUtils_hh
#define G4InclUtils_hh 1

#include "globals.hh"

class G4InclUtils
{
 public:
  static G4double calculate4MomentumScaling(G4int A, G4int Z, G4double excitationE, G4double kineticE,
					    G4double px, G4double py, G4double pz);

protected:
  G4InclUtils();
  ~G4InclUtils();

public:
  // Verbosity levels:
  static const G4int silent = 0;
  static const G4int debug = 2;
  static const G4int verbose = 3;
  static const G4int verboseAll = 5;
};

#endif
