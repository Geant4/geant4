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
// $Id: G4ionIonisation.hh 106717 2017-10-20 09:41:27Z gcosmo $
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
// 12-09-08 Removed InitialiseMassCharge and CorrectionsAlongStep (VI)
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

class G4Material;
class G4EmCorrections;

class G4ionIonisation : public G4VEnergyLossProcess
{
public:

  explicit G4ionIonisation(const G4String& name = "ionIoni");

  virtual ~G4ionIonisation();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) final;

  void AddStoppingData(G4int Z, G4int A, const G4String& materialName,
		       G4PhysicsVector* dVector);

  void ActivateStoppingData(G4bool);

  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  // Print out of the class parameters
  virtual void StreamProcessInfo(std::ostream& outFile,
                             G4String endOfLine=G4String("\n")) const override;

  virtual void 
  InitialiseEnergyLossProcess(const G4ParticleDefinition*,
			      const G4ParticleDefinition*) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, G4double cut) final;

  inline G4double BetheBlochEnergyThreshold();

private:

  // hide assignment operator
  G4ionIonisation & operator=(const G4ionIonisation &right) = delete;
  G4ionIonisation(const G4ionIonisation&) = delete;

  G4EmCorrections*            corr;

  const G4ParticleDefinition* theParticle;

  G4double   eth;

  G4bool     isInitialised;
  G4bool     stopDataActive;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ionIonisation::ActivateStoppingData(G4bool val)
{
  stopDataActive = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionIonisation::BetheBlochEnergyThreshold()
{
  return eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
