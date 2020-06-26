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
// G4PhysicsModelCatalog
//
// Class description:
//
// Singleton, collection of physics models, to be used by models and G4Track.

// Author: M.Asai (SLAC), 26 September 2013
// --------------------------------------------------------------------
#ifndef G4PhysicsModelCatalog_hh
#define G4PhysicsModelCatalog_hh

#include <vector>

#include "G4String.hh"
#include "globals.hh"

class G4PhysicsModelCatalog
{
 public:
  ~G4PhysicsModelCatalog();
  G4PhysicsModelCatalog(const G4PhysicsModelCatalog&) = delete;
  G4PhysicsModelCatalog& operator=(const G4PhysicsModelCatalog&) = delete;

  static G4int Register(const G4String&);
  static const G4String& GetModelName(G4int);

  static G4int GetIndex(const G4String&);
  static G4int Entries();
  static void Destroy();

 private:
  G4PhysicsModelCatalog();

  static std::vector<G4String>* theCatalog;
};

#endif
