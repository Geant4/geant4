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
// $Id: G4ionIonisation.hh,v 1.50 2007/11/09 11:45:45 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ionIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 26-12-02 Secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 15-02-03 Add control on delta pointer (V.Ivanchenko)
// 23-05-03 Add fluctuation model as a member function (V.Ivanchenko)
// 03-08-03 Add effective charge and saturation of tmax (V.Ivanchenko)
// 12-11-03 Fix problem of negative effective charge (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 11-04-05 Move MaxSecondary energy to model (V.Ivanchneko)
// 11-04-04 Move MaxSecondaryEnergy to models (V.Ivanchenko)
// 10-05-06 Add a possibility to download user data (V.Ivantchenko)
// 22-07-06 Remove obsolete method (V.Ivantchenko)
// 07-11-07 Moved CorrectionsAlongStep to cc (V.Ivantchenko)
//
// Class Description:
//
// This class manages the ionisation process for ions. Effective charge,
// nuclear stopping power, energy loss corrections are taken into account.
// It inherites from G4VEnergyLossLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4ionIonisation_h
#define G4ionIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4VEmModel.hh"
#include "G4EmCorrections.hh"

class G4Material;
class G4PhysicsVector;
class G4BraggIonModel;

class G4ionIonisation : public G4VEnergyLossProcess
{
public:

  G4ionIonisation(const G4String& name = "ionIoni");

  virtual ~G4ionIonisation();

  inline G4bool IsApplicable(const G4ParticleDefinition& p);

  // Print out of the class parameters
  virtual void PrintInfo();

  void AddStoppingData(G4int Z, G4int A, const G4String& materialName,
		       G4PhysicsVector& dVector);

  void ActivateStoppingData(G4bool);

  void ActivateNuclearStopping(G4bool);

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

  virtual void CorrectionsAlongStep(const G4MaterialCutsCouple*,
				    const G4DynamicParticle*,
				    G4double& eloss,
				    G4double& length);

  inline void InitialiseMassCharge(const G4Track&);

  inline G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				   const G4Material*, G4double cut);

  inline G4double BetheBlochEnergyThreshold();

  inline G4bool NuclearStoppingFlag();

  // protected pointers 
  G4ionEffectiveCharge*       effCharge;
  G4EmCorrections*            corr;

private:

  // hide assignment operator
  G4ionIonisation & operator=(const G4ionIonisation &right);
  G4ionIonisation(const G4ionIonisation&);

  // cash
  const G4Material*           curMaterial;
  const G4ParticleDefinition* curParticle;
  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

  G4double                    preKinEnergy;

  G4double                    eth;
  G4double                    baseMass;
  G4double                    massRatio;
  G4double                    massFactor;
  G4double                    charge2;

  G4bool                      isInitialised;
  G4bool                      stopDataActive;
  G4bool                      nuclearStopping;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4ionIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived() &&
          p.GetParticleType() == "nucleus");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionIonisation::MinPrimaryEnergy(
          const G4ParticleDefinition*, const G4Material*, G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double g = std::sqrt(1. + x);
  return proton_mass_c2*(g - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ionIonisation::InitialiseMassCharge(const G4Track& track)
{
  preKinEnergy = track.GetKineticEnergy();
  massRatio    = baseMass/track.GetDynamicParticle()->GetMass();
  charge2      = effCharge->EffectiveChargeSquareRatio(track.GetDefinition(),
						       track.GetMaterial(),
						       preKinEnergy);
  SetDynamicMassCharge(massRatio, charge2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ionIonisation::ActivateStoppingData(G4bool val)
{
  stopDataActive = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ionIonisation::ActivateNuclearStopping(G4bool val)
{
  nuclearStopping = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionIonisation::BetheBlochEnergyThreshold()
{
  return eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4ionIonisation::NuclearStoppingFlag()
{
  return nuclearStopping;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
