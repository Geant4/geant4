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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for registration of static cross sections components
//
// Author: V.Ivanchenko 31.05.2018
//
// Modifications:
// 07.03.2024 V.Ivanchenko updated signature - now it is a store without any access
//                         to stored objects
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronXSDataTable_h
#define G4HadronXSDataTable_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4PiData.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4HadronXSDataTable
{

public:

  static G4HadronXSDataTable* Instance();

  ~G4HadronXSDataTable();

  void AddPiData(std::vector<G4PiData*>* ptr);

  void AddTable(G4PhysicsTable* ptr);

  // Assignment operator and copy constructor
  G4HadronXSDataTable & operator=(const G4HadronXSDataTable &right) = delete;
  G4HadronXSDataTable(const G4HadronXSDataTable&) = delete;

private:

  G4HadronXSDataTable();

  static G4HadronXSDataTable* sInstance;

  std::vector<std::vector<G4PiData*>* > fPiData;
  std::vector<G4PhysicsTable*> fTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
