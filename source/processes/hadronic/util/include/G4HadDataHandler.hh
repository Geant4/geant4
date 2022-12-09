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
// File name:     G4HadDataHandler
//
// Author:        V. Ivanchenko 
// 
// Creation date: 05 July 2022
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

#ifndef G4HadDataHandler_h
#define G4HadDataHandler_h 1

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VProcess;

class G4HadDataHandler
{
public:

  explicit G4HadDataHandler(std::size_t nTable);

  ~G4HadDataHandler();

  // add table
  void AddTable(G4PhysicsTable*);

  // update table with existing index
  void UpdateTable(G4PhysicsTable*, std::size_t idx);

  // clean existing table
  void CleanTable(std::size_t idx);
  
  // keep pointer of master thread process
  void SetMasterProcess(const G4VProcess*);

  // return pointer of master thread process
  const G4VProcess* GetMasterProcess(std::size_t idx) const;

  // access to the table via index
  inline const G4PhysicsTable* GetTable(std::size_t idx) const { 
    return (idx < tLength) ? data[idx] : nullptr; 
  }
  inline G4PhysicsTable* Table(std::size_t idx) const {
    return (idx < tLength) ? data[idx] : nullptr; 
  }

  // access to vector, no check on input index
  inline 
  const G4PhysicsVector* GetVector(std::size_t itable, std::size_t ivec) const
  { return (*(data[itable]))[ivec]; }

  inline const std::vector<G4PhysicsTable*>& GetTables() const { return data; }

  // hide assignment operator 
  G4HadDataHandler & operator=(const G4HadDataHandler &right) = delete;
  G4HadDataHandler(const G4HadDataHandler&) = delete;

private:

  std::vector<G4PhysicsTable*> data;
  std::size_t tLength;
  std::vector<const G4VProcess*> masterProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

