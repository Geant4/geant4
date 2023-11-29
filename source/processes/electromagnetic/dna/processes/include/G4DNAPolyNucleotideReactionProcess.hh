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
/// file: G4DNAPolyNucleotideReactionProcess.hh
/// brief: This file handls reaction process with DNA geometry.
#ifndef G4DNAPolyNucleotideReactionProcess_hh
#define G4DNAPolyNucleotideReactionProcess_hh
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4VITDiscreteProcess.hh"
#include <variant>
class G4DNAMolecularReactionTable;
class G4DNAComponentNode;
class G4VDNAHitModel;

class G4DNAPolyNucleotideReactionProcess : public G4VITDiscreteProcess
{
 public:
  explicit G4DNAPolyNucleotideReactionProcess(
    const G4String& aName = "DNAStaticMoleculeReactionProcess",
    G4int verbosityLevel  = 0);
  ~G4DNAPolyNucleotideReactionProcess() override;

  inline void SetDNADamageReactionModel(G4VDNAHitModel* pModel);

  G4bool IsApplicable(const G4ParticleDefinition&) override { return true; }

  G4double CalculateTimeStep(const G4Track& trackA,
                             const G4double& userTimeStep = 0);

  void StartTracking(G4Track* aTrack) override;

  G4double PostStepGetPhysicalInteractionLength(
    const G4Track&,  // track
    G4double,        // previousStepSize
    G4ForceCondition* pForceCond) override;
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step&) override;

  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override
  {
    return DBL_MAX;
  }

  G4DNAPolyNucleotideReactionProcess& operator =(
    const G4DNAPolyNucleotideReactionProcess&) = delete;
  inline void SetVerbose(G4int verbose);

 protected:
  G4VParticleChange fParticleChange;

 private:
  struct G4PolyNucleotideReactionState
    : public G4ProcessStateBase<G4PolyNucleotideReactionState>
  {
    G4PolyNucleotideReactionState();
    ~G4PolyNucleotideReactionState() override = default;
    G4String GetType() override { return "PolyNucleotideReactionState"; }

    using DNANode =
      std::variant<const G4DNAComponentNode*, /*for dnadamage chain*/
                   const G4VPhysicalVolume* /*for molecularDNA chain*/>;
    DNANode fNodeReactant;
    G4double fSampledMinTimeStep;
    G4double fPreviousTimeAtPreStepPoint;
  };
  G4bool fHasAlreadyReachedNullTime;
  G4int fVerbose;
  G4double fRCutOff;
  G4VDNAHitModel* fpDamageModel;
  G4DNAPolyNucleotideReactionProcess(const G4DNAPolyNucleotideReactionProcess&);
};
inline void G4DNAPolyNucleotideReactionProcess::SetVerbose(G4int verbose)
{
  fVerbose = verbose;
}

inline void G4DNAPolyNucleotideReactionProcess::SetDNADamageReactionModel(
  G4VDNAHitModel* pModel)
{
  fpDamageModel = pModel;
}
#endif
