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
// $Id: G4EmDataHandler.hh 73844 2013-09-13 14:16:30Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4EmDataHandler
//
// Author:        V. Ivanchenko 
// 
// Creation date: 16 August 2016
//
// Modifications: 
//
// Class Description: 
//
// A storage of G4PhysicsTable
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4EmDataHandler_h
#define G4EmDataHandler_h 1

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;

class G4EmDataHandler
{
public:

  explicit G4EmDataHandler(size_t nTable);

  ~G4EmDataHandler();

  size_t SetTable(G4PhysicsTable*);

  G4PhysicsTable* MakeTable();

  G4bool StorePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& namet, 
			   const G4String& procname, 
			   G4bool ascii);

  G4bool RetrievePhysicsTable(size_t idx,
			      const G4ParticleDefinition* part,
			      const G4String& namet, 
			      const G4String& procname, 
			      G4bool ascii);

  inline G4PhysicsTable* GetTable(size_t idx) { return data[idx]; }

  inline std::vector<G4PhysicsTable*>& GetTables() { return data; }

  // hide assignment operator 
  G4EmDataHandler & operator=(const G4EmDataHandler &right) = delete;
  G4EmDataHandler(const G4EmDataHandler&) = delete;

private:

  std::vector<G4PhysicsTable*> data;
  size_t tLength;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

