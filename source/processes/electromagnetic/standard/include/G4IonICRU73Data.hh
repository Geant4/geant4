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

#ifndef G4IonICRU73Data_h
#define G4IonICRU73Data_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4IonICRU73Data
//
// Description: Data on stopping power
//
// Author:      Vladimir Ivanchenko
//
// Creation date: 23.10.2021
//
// Class Description:
//
// Container for parameterised data on ion stopping power
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <vector>
#include <iostream>
#include "globals.hh"

class G4PhysicsLogVector;
class G4PhysicsFreeVector;
class G4Material;

class G4IonICRU73Data 
{ 
public: 

  explicit G4IonICRU73Data();

  ~G4IonICRU73Data();

  // should be called before each run
  void Initialise();

  // Z - atomic number of the projectile
  // e - energy per nucleon, loge = log(e)
  G4double GetDEDX(const G4Material*, const G4int Z,
                   const G4double e, const G4double loge) const;

  // hide assignment operator
  G4IonICRU73Data & operator = (const  G4IonICRU73Data &right) = delete;
  G4IonICRU73Data(const  G4IonICRU73Data&) = delete;

private:

  void ReadMaterialData(const G4Material* mat, const G4double fact, 
                        const G4bool type);

  void ReadElementData(const G4Material* mat, G4bool type);

  G4PhysicsLogVector* FindOrBuildElementData(const G4int Z, 
                                             const G4int Z1, 
                                             G4bool useICRU90);

  G4PhysicsLogVector* RetrieveVector(std::ostringstream& in,
                                     G4bool warn);

  G4double fEmin;
  G4double fEmax;

  std::vector<G4int> fMatIndex;
  // projectile Z <= 80, target element Z <= 92
  std::vector<G4PhysicsLogVector*>* fMatData[81] = {nullptr};
  G4PhysicsLogVector* fElmData[81][93] = {{nullptr}};
  G4PhysicsFreeVector* fVector = nullptr;

  G4int fNmat = 0;
  G4int fNbins = 0;
  G4int fNbinsPerDecade = 10;
  G4int fVerbose = 0;
  G4bool fSpline = false;
  G4String fDataDirectory = "";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
