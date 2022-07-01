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
#include "G4PhysicsVector.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmProcess;
class G4VEnergyLossProcess;

class G4EmDataHandler
{
public:

  explicit G4EmDataHandler(size_t nTable);

  ~G4EmDataHandler();

  // add table
  size_t SetTable(G4PhysicsTable*);

  // update existing table
  void UpdateTable(G4PhysicsTable*, size_t idx);

  // assuming that the table is already defined
  G4PhysicsTable* MakeTable(size_t idx);

  // existing table may be substituted
  G4PhysicsTable* MakeTable(G4PhysicsTable*, size_t idx);

  // clean existing table 
  void CleanTable(size_t idx);

  G4bool StorePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& fname, 
			   G4bool ascii);

  G4bool RetrievePhysicsTable(size_t idx,
			      const G4ParticleDefinition* part,
			      const G4String& fname,
			      G4bool ascii, G4bool spline);
  
  void SetMasterProcess(const G4VEmProcess*);

  const G4VEmProcess* GetMasterProcess(size_t idx) const;

  inline const G4PhysicsTable* GetTable(size_t idx) const { 
    return (idx < tLength) ? data[idx] : nullptr; 
  }

  inline G4PhysicsTable* Table(size_t idx) const {
    return (idx < tLength) ? data[idx] : nullptr; 
  }

  inline const G4PhysicsVector* GetVector(size_t itable, size_t ivec) const
  { return (*(data[itable]))[ivec]; }

  inline const std::vector<G4PhysicsTable*>& GetTables() const { return data; }

  // hide assignment operator 
  G4EmDataHandler & operator=(const G4EmDataHandler &right) = delete;
  G4EmDataHandler(const G4EmDataHandler&) = delete;

private:

  std::vector<G4PhysicsTable*> data;
  size_t tLength;
  std::vector<const G4VEmProcess*> masterProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

