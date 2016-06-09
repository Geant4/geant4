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
// $Id: G4ionGasIonisation.hh,v 1.2 2007/08/13 06:13:30 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ionGasIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.07.2007
//
// Modifications:
//
//
// Class Description:
//
// This class manages the ionisation process for ions. It inherites 
// from G4ionIonisation and sample ion/media charge exchange
//

// -------------------------------------------------------------------
//

#ifndef G4ionGasIonisation_h
#define G4ionGasIonisation_h 1

#include "G4ionIonisation.hh"

class G4Material;

class G4ionGasIonisation : public G4ionIonisation
{
public:

  G4ionGasIonisation(const G4String& name = "ionGasIoni");

  virtual ~G4ionGasIonisation();

  virtual void PrintInfo();

protected:

  // Initialise process before run
  virtual void InitialiseEnergyLossProcess(
				   const G4ParticleDefinition*,
                                   const G4ParticleDefinition*);

  // Initialise dynamic charge before step
  virtual void InitialiseMassCharge(const G4Track&);

  // Apply correction after step and modify dynamic charge
  virtual void CorrectionsAlongStep(
                           const G4MaterialCutsCouple*,
                           const G4DynamicParticle*,
                                 G4double& eloss,
                                 G4double& length);

private:

  // Sample change of charge of the projectile ion
  G4double SampleChargeAfterStep(G4double qeff, G4double xeff);

  // hide assignment operator
  G4ionGasIonisation & operator=(const G4ionGasIonisation &right);
  G4ionGasIonisation(const G4ionGasIonisation&);

  const G4ParticleDefinition* currParticle;
  const G4ParticleDefinition* baseParticle;

  G4double                    basePartMass;
  G4double                    currMassRatio;
  G4double                    atomXS;
  G4double                    preStepKinEnergy;

  G4int                       currentIonZ;
  G4int                       ionZ;

  G4bool                      initialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
