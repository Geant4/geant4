//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4UniversalFluctuation
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.01.2002
//
// Modifications: 09.12.2002 VI remove warnings
//
// Class Description: 
//
// Implementation of energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuation_h
#define G4UniversalFluctuation_h 1


#include "G4VEmFluctuationModel.hh"

class G4UniversalFluctuation : public G4VEmFluctuationModel
{

public:

  G4UniversalFluctuation();

  ~G4UniversalFluctuation();

  virtual void SampleFluctuations(const G4Material*, const G4DynamicParticle*,
 				        G4double&, G4double&, G4double); 

  virtual void Initialise(const G4ParticleDefinition*);

protected:

private:

  // hide assignment operator 
  G4UniversalFluctuation & operator=(const  G4UniversalFluctuation &right);
  G4UniversalFluctuation(const  G4UniversalFluctuation&);

  const G4ParticleDefinition* particle;
  G4double particleMass;
  G4double chargeSquare;

  const G4Material* lastMaterial;
  G4int materialIndex;

  // data members to speed up the fluctuation calculation
  G4double ipotFluct;
  G4double electronDensity;
  G4double zeff;
  
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double rateFluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double ipotLogFluct;

  G4double minNumberInteractionsBohr;
  G4double theBohrBeta2;
  G4double minLoss;
  G4double problim;
  G4double sumalim;
  G4double alim;
  G4int    nmaxCont1;
  G4int    nmaxCont2;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

