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
// File name:     G4PolarizedCompton
//
// Author:        Andreas Schaelicke
//                based on code by Michel Maire / Vladimir IVANTCHENKO
//
// Class description
//   polarized version of Compton scattering
//
// -----------------------------------------------------------------------------

#ifndef PolarizedComptonScattering_h
#define PolarizedComptonScattering_h 1

#include "globals.hh"
#include "G4Gamma.hh"
#include "G4VEmProcess.hh"

class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4PolarizedComptonModel;

class G4PolarizedCompton : public G4VEmProcess

{
 public:
  explicit G4PolarizedCompton(const G4String& processName = "pol-compt",
                              G4ProcessType type          = fElectromagnetic);

  virtual ~G4PolarizedCompton() override;

  // true for Gamma only.
  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;

  virtual void ProcessDescription(std::ostream&) const override;
  virtual void DumpInfo() const override { ProcessDescription(G4cout); };

  void SetModel(const G4String& name);

  G4PolarizedCompton& operator=(const G4PolarizedCompton& right) = delete;
  G4PolarizedCompton(const G4PolarizedCompton&)                  = delete;

 protected:
  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  virtual G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

 private:
  static G4PhysicsTable* theAsymmetryTable;  // table for crosssection asymmetry
  void CleanTable();

  void BuildAsymmetryTable(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy, const G4MaterialCutsCouple* couple,
                            const G4ParticleDefinition& particle, G4double cut,
                            G4double& tAsymmetry);

  G4double ComputeSaturationFactor(const G4Track& aTrack);

  G4PolarizedComptonModel* fEmModel;

  G4int fType;

  G4bool fBuildAsymmetryTable;
  G4bool fUseAsymmetryTable;

  G4bool fIsInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
