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
// File name:     G4IonFluctuation
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 02.04.2003
//
// Modifications:
//
// Class Description:
//
// Implementation of ion energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4IonFluctuation_h
#define G4IonFluctuation_h 1


#include "G4VEmFluctuationModel.hh"

class G4IonFluctuation : public G4VEmFluctuationModel
{

public:

  G4IonFluctuation(const G4String& nam = "IonFluc");

  ~G4IonFluctuation();

  G4double SampleFluctuations(const G4Material*,
                              const G4DynamicParticle*,
                                    G4double&,
                                    G4double&,
                                    G4double&);

  G4double Dispersion(    const G4Material*,
                          const G4DynamicParticle*,
 				G4double&,
                                G4double&);

  void Initialise(const G4ParticleDefinition*);

private:

  G4double CoeffitientA(G4double&);
  G4double CoeffitientB(const G4Material*, G4double&);

  // hide assignment operator
  G4IonFluctuation & operator=(const  G4IonFluctuation &right);
  G4IonFluctuation(const  G4IonFluctuation&);

  const G4ParticleDefinition* particle;

  G4double particleMass;
  G4double charge;
  G4double chargeSquare;

  // data members to speed up the fluctuation calculation
  G4double minNumberInteractionsBohr;
  G4double theBohrBeta2;
  G4double minFraction;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

