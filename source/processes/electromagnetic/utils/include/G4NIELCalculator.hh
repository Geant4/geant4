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
#ifndef G4NIELCalculator_h
#define G4NIELCalculator_h 1

// -------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4NIELCalculator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 30.05.2019
//
// Modifications:
//
//
// Class Description:
// This is a helper class to compute NIEL in user stepping action 
// or sensitive detector code. User should provide G4VEmModel
// objects, which has NIEL model
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4VEmModel;

class G4NIELCalculator
{
public: 

  explicit G4NIELCalculator(G4VEmModel*, G4int verb);

  ~G4NIELCalculator() = default;

  // initialisation before start of run
  void Initialise();

  // compute nuclear stopping power
  G4double ComputeNIEL(const G4Step*);

  // kinetic energy of recoil nucleus or zero
  G4double RecoilEnergy(const G4Step*);

  // replace model of NIEL
  void AddEmModel(G4VEmModel*); 

  // hide assignment operator
  G4NIELCalculator & operator=(const G4NIELCalculator &right) = delete;
  G4NIELCalculator(const G4NIELCalculator&) = delete;

private:

  G4VEmModel* fModel;  
  G4int fVerbose;             
};

#endif

