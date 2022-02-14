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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4ContinuousGainOfEnergy
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Continuous process acting on adjoint particles to compute the continuous
//  gain of energy of charged particles when they are tracked back.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4ContinuousGainOfEnergy_h
#define G4ContinuousGainOfEnergy_h 1

#include "globals.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VContinuousProcess.hh"

class G4Material;
class G4MaterialCutsCouple;
class G4ParticleChange;
class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VEmModel;
class G4VEnergyLossProcess;

class G4ContinuousGainOfEnergy : public G4VContinuousProcess
{
 public:
  explicit G4ContinuousGainOfEnergy(const G4String& name = "EnergyGain",
                                    G4ProcessType type   = fElectromagnetic);

  ~G4ContinuousGainOfEnergy() override;

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  void SetLossFluctuations(G4bool val);

  inline void SetDirectEnergyLossProcess(G4VEnergyLossProcess* aProcess)
  {
    fDirectEnergyLossProcess = aProcess;
  };

  void SetDirectParticle(G4ParticleDefinition* p);

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };

  G4ContinuousGainOfEnergy(G4ContinuousGainOfEnergy&) = delete;
  G4ContinuousGainOfEnergy& operator=(const G4ContinuousGainOfEnergy& right) =
    delete;

 protected:
  G4double GetContinuousStepLimit(const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimumStep,
                                  G4double& currentSafety) override;

 private:
  void DefineMaterial(const G4MaterialCutsCouple* couple);
  void SetDynamicMassCharge(const G4Track& track, G4double energy);

  const G4Material* fCurrentMaterial = nullptr;
  const G4MaterialCutsCouple* fCurrentCouple = nullptr;

  G4VEmModel* fCurrentModel                      = nullptr;
  G4VEnergyLossProcess* fDirectEnergyLossProcess = nullptr;
  G4ParticleDefinition* fDirectPartDef           = nullptr;

  G4double fCurrentTcut      = 0.;
  G4double fPreStepKinEnergy = 1.;
  G4double fLinLossLimit     = 0.05;
  G4double fMassRatio        = 1.;

  size_t fCurrentCoupleIndex = 9999999;

  G4bool fIsIon                      = false;
  G4bool fLossFluctuationFlag        = true;
  G4bool fLossFluctuationArePossible = true;
};

///////////////////////////////////////////////////////
inline void G4ContinuousGainOfEnergy::DefineMaterial(
  const G4MaterialCutsCouple* couple)
{
  if(couple != fCurrentCouple)
  {
    fCurrentCouple      = couple;
    fCurrentMaterial    = couple->GetMaterial();
    fCurrentCoupleIndex = couple->GetIndex();

    const std::vector<G4double>* aVec =
      G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
    fCurrentTcut = (*aVec)[fCurrentCoupleIndex];
  }
}

#endif
