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
// $Id: G4EmSaturation.hh 100346 2016-10-18 15:30:36Z gcosmo $
//
//
#ifndef G4EmSaturation_h
#define G4EmSaturation_h 1

// -------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmSaturation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.02.2008
//
// Modifications:
//
//
// Class Description:
//   Compution on saturation effect, which reduce visible energy 
//   deposition at the step. Default implementation takes into 
//   account Birks effect. Birks coefficients for some materials
//   from G4 database on materials are provided
//
//   This class assumed to be G4ThreadLocal, because it is using 
//   cache value for material 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4NistManager;
class G4MaterialCutsCouple;
class G4Material;

class G4EmSaturation
{
public: 

  explicit G4EmSaturation(G4int verb);
  virtual ~G4EmSaturation();

  // this method may be overwritten in the derived class
  // which implements alternative algorithm of saturation
  virtual G4double VisibleEnergyDeposition(const G4ParticleDefinition*, 
					   const G4MaterialCutsCouple*,
					   G4double length, 
					   G4double edepTotal,
					   G4double edepNIEL = 0.0) const;

  // activate default model  
  void InitialiseG4Saturation();

  // find and Birks coefficient 
  G4double FindG4BirksCoefficient(const G4Material*);

  // dump coeffitients used in run time
  void DumpBirksCoefficients();

  // dump G4 list
  void DumpG4BirksCoefficients();

  // this method should not be overwitten
  inline G4double VisibleEnergyDepositionAtAStep(const G4Step*) const; 

  inline void SetVerbose(G4int);

private:

  // hide assignment operator
  G4EmSaturation & operator=(const G4EmSaturation &right) = delete;
  G4EmSaturation(const G4EmSaturation&) = delete;

  void InitialiseBirksCoefficient(const G4Material*);

  void InitialiseG4materials();

  const G4ParticleDefinition* electron;
  const G4ParticleDefinition* proton;
  G4NistManager*              nist;

  G4int    verbose;             
  G4int    nG4Birks;
  G4int    nWarnings;

  static G4int nMaterials;

  // list of materials used in run time
  static std::vector<G4double> massFactors;
  static std::vector<G4double> effCharges;

  // list of G4 materials 
  static std::vector<G4double>  g4MatData;
  static std::vector<G4String>  g4MatNames;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EmSaturation::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EmSaturation::VisibleEnergyDepositionAtAStep(
                const G4Step* step) const
{
  return VisibleEnergyDeposition(step->GetTrack()->GetParticleDefinition(),
                                 step->GetTrack()->GetMaterialCutsCouple(),
				 step->GetStepLength(),
                                 step->GetTotalEnergyDeposit(),
                                 step->GetNonIonizingEnergyDeposit());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

