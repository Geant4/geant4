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
// $Id: G4MuIonisation.hh 106716 2017-10-20 09:40:09Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MuIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 30.05.1997
//
// Modifications:
//
// corrected by L.Urban on 24/09/97
// corrected by L.Urban on 13/01/98
// bugs fixed by L.Urban on 02/02/99
// 10/02/00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 19-09-01 come back to previous process name "hIoni"
// 29-10-01 all static functions no more inlined
// 10-05-02 V.Ivanchenko update to new design
// 09-12-02 V.Ivanchenko remove warning
// 26-12-02 Secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Add fluctuation model as a member function (V.Ivanchenko)
// 03-06-03 Add SetIntegral method to choose fluctuation model (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 17-08-04 Rename the process "Mu" -> "mu" (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//
// Class Description:
//
// This class manages the ionisation process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossProcess.
//

// -------------------------------------------------------------------
//

#ifndef G4MuIonisation_h
#define G4MuIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"
#include "G4VEmModel.hh"

class G4Material;

class G4MuIonisation : public G4VEnergyLossProcess
{

public:

  explicit G4MuIonisation(const G4String& name = "muIoni");

  virtual ~G4MuIonisation();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, G4double cut) override;

  // Print out of the class parameters
  virtual void PrintInfo() override;

  // print description in html
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                               const G4ParticleDefinition*) override;

private:

  // hide assignment operator
  G4MuIonisation & operator=(const G4MuIonisation &right) = delete;
  G4MuIonisation(const G4MuIonisation&) = delete;

  G4double    mass;
  G4double    ratio;

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;
  G4bool                      isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
