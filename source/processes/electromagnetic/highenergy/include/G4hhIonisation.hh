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
// $Id: G4hhIonisation.hh,v 1.3 2006/10/27 14:11:12 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4hhIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 30.09.2005
//
// Modifications:
//
//
// Class Description:
//
// This class manages the ionisation process for exsotic hadrons
// without any delta electrons production
//

// -------------------------------------------------------------------
//

#ifndef G4hhIonisation_h
#define G4hhIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "globals.hh"
#include "G4VEmModel.hh"

class G4Material;
class G4VEmFluctuationModel;

class G4hhIonisation : public G4VEnergyLossProcess
{

public:

  G4hhIonisation(const G4String& name = "hhIoni");

  virtual ~G4hhIonisation();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
			    const G4Material*, G4double cut);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  std::vector<G4DynamicParticle*>*  SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&);

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

private:

  // hide assignment operator
  G4hhIonisation & operator=(const G4hhIonisation &right);
  G4hhIonisation(const G4hhIonisation&);

  G4double   mass;
  G4double   ratio;
  G4double   minKinEnergy;
  G4double   maxKinEnergy;

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;
  G4VEmFluctuationModel*      flucModel;

  G4bool                      isInitialised;

  G4double                    eth;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4hhIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 100.0*MeV &&
	 !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hhIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
						const G4Material*,
						G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double y = electron_mass_c2/mass;
  G4double g = x*y + std::sqrt((1. + x)*(1. + x*y*y));
  return mass*(g - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4hhIonisation::SecondariesPostStep(
                                                  G4VEmModel*,
                                            const G4MaterialCutsCouple*,
                                            const G4DynamicParticle*,
                                                  G4double&)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
