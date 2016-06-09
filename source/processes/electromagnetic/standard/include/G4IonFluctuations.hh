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
// $Id: G4IonFluctuations.hh,v 1.1 2004/12/01 17:36:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4IonFluctuations
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications: 
//
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
//
// Class Description:
//
// Implementation of ion energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4IonFluctuations_h
#define G4IonFluctuations_h 1


#include "G4VEmFluctuationModel.hh"

class G4IonFluctuations : public G4VEmFluctuationModel
{

public:

  G4IonFluctuations(const G4String& nam = "IonFluc");

  virtual ~G4IonFluctuations();

  G4double SampleFluctuations(const G4Material*,
                              const G4DynamicParticle*,
                                    G4double&,
                                    G4double&,
                                    G4double&);

  G4double Dispersion(    const G4Material*,
                          const G4DynamicParticle*,
 				G4double&,
                                G4double&);

  void InitialiseMe(const G4ParticleDefinition*);

private:

  G4double CoeffitientA(G4double&);
  G4double CoeffitientB(const G4Material*, G4double&);
  G4double RelativisticFactor(const G4Material*, G4double&);

  // hide assignment operator
  G4IonFluctuations & operator=(const  G4IonFluctuations &right);
  G4IonFluctuations(const  G4IonFluctuations&);

  const G4ParticleDefinition* particle;

  G4double particleMass;
  G4double charge;
  G4double chargeSquare;
  G4double chargeSqRatio;

  // data members to speed up the fluctuation calculation
  G4double minNumberInteractionsBohr;
  G4double theBohrBeta2;
  G4double minFraction;
  G4double xmin;
  G4double minLoss;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

