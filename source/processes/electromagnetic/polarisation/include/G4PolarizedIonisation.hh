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
// -------------------------------------------------------------------
//
// Geant4 Class header file
//
// File name:     G4PolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Class Description:
//   polarized version of G4eIonisation
// ----------------------------------------------------------------------------

#ifndef G4PolarizedIonisation_h
#define G4PolarizedIonisation_h 1

#include "globals.hh"
#include "G4VEnergyLossProcess.hh"

class G4Material;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4VEmFluctuationModel;
class G4PolarizedIonisationModel;
class G4Track;

class G4PolarizedIonisation : public G4VEnergyLossProcess
{
 public:
  explicit G4PolarizedIonisation(const G4String& name = "pol-eIoni");

  virtual ~G4PolarizedIonisation() override;

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override;

  virtual void ProcessDescription(std::ostream&) const override;
  virtual void DumpInfo() const override { ProcessDescription(G4cout); };

  G4PolarizedIonisation& operator=(const G4PolarizedIonisation& right) = delete;
  G4PolarizedIonisation(const G4PolarizedIonisation&)                  = delete;

 protected:
  virtual void InitialiseEnergyLossProcess(
    const G4ParticleDefinition*, const G4ParticleDefinition*) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut) override;

  virtual G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

  virtual G4double GetMeanFreePath(const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

 private:
  void CleanTables();

  void BuildAsymmetryTables(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy, const G4MaterialCutsCouple* couple,
                            const G4ParticleDefinition& particle, G4double cut,
                            G4double& tasm);

  G4double ComputeSaturationFactor(const G4Track& aTrack);

  G4VEmFluctuationModel* fFlucModel;
  G4PolarizedIonisationModel* fEmModel;

  G4PhysicsTable* fAsymmetryTable;
  G4PhysicsTable* fTransverseAsymmetryTable;

  G4bool fIsElectron;
  G4bool fIsInitialised;
};

#endif
