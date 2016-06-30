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
// $Id: G4ePolarizedIonisation.hh 96114 2016-03-16 18:51:33Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ePolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 10-11-05, include polarization description (A.Schaelicke)
// , create asymmetry table and determine interactionlength 
// , update polarized differential cross section 
//
// 20-08-06, modified interface (A.Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
// Class Description:
//
// polarized version of G4eIonisation
// ----------------------------------------------------------------------------
// -------------------------------------------------------------------
//

#ifndef G4ePolarizedIonisation_h
#define G4ePolarizedIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;
class G4ParticleDefinition;
class G4VEmFluctuationModel;
class G4PolarizedMollerBhabhaModel;

class G4ePolarizedIonisation : public G4VEnergyLossProcess
{

public:

  explicit G4ePolarizedIonisation(const G4String& name = "pol-eIoni");

  virtual ~G4ePolarizedIonisation();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override;

  // Print out of the class parameters
  virtual void PrintInfo() override;

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut) override;

  // for polarization
  virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) override;

  virtual G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

private:

  void CleanTables();

  void BuildAsymmetryTables(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& particle,
			    G4double cut,
			    G4double &tasm);

  G4double ComputeSaturationFactor(const G4Track& aTrack);

  G4ePolarizedIonisation & 
    operator=(const G4ePolarizedIonisation &right) = delete;
  G4ePolarizedIonisation(const G4ePolarizedIonisation&) = delete;

  G4ParticleDefinition* theElectron;
  G4VEmFluctuationModel* flucModel;
  G4PolarizedMollerBhabhaModel* emModel;

  G4bool isElectron;
  G4bool isInitialised;

  // for polarization:
  G4ThreeVector theTargetPolarization;
  
  G4PhysicsTable* theAsymmetryTable;  
  G4PhysicsTable* theTransverseAsymmetryTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
