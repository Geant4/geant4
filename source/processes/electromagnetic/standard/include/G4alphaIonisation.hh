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
// $Id: G4alphaIonisation.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4alphaIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.10.2009
//
// Modifications:
//
//
// Class Description:
//
// This class manages the ionisation process for He ions. Effective charge,
// nuclear stopping power, energy loss corrections are taken into account.
// It inherites from G4VEnergyLossLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4alphaIonisation_h
#define G4alphaIonisation_h 1

#include "G4VEnergyLossProcess.hh"

class G4Material;

class G4alphaIonisation : public G4VEnergyLossProcess
{
public:

  explicit G4alphaIonisation(const G4String& name = "alphaIoni");

  virtual ~G4alphaIonisation();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) final;

  // Print out of the class parameters
  virtual void PrintInfo() override;

  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void 
  InitialiseEnergyLossProcess(const G4ParticleDefinition*,
			      const G4ParticleDefinition*) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				   const G4Material*, G4double cut) final;

  inline G4double BetheBlochEnergyThreshold();

private:

  // hide assignment operator
  G4alphaIonisation & operator=(const G4alphaIonisation &right) = delete;
  G4alphaIonisation(const G4alphaIonisation&) = delete;

  const G4ParticleDefinition* theParticle;

  G4double   mass;
  G4double   ratio;
  G4double   eth;

  G4bool     isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4alphaIonisation::BetheBlochEnergyThreshold()
{
  return eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
