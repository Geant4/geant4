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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4NeutronGeneralProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 08.08.2022
//
// Modifications:
//
// Class Description:
//
// It is the neutron super process

// -------------------------------------------------------------------
//

#ifndef G4NeutronGeneralProcess_h
#define G4NeutronGeneralProcess_h 1

#include "G4HadronicProcess.hh"
#include "globals.hh"
#include "G4HadDataHandler.hh"
#include <vector>

class G4Step;
class G4Track;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VCrossSectionDataSet;
class G4CrossSectionDataStore;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4NeutronGeneralProcess : public G4HadronicProcess
{
public:

  explicit G4NeutronGeneralProcess(const G4String& pname="NeutronGeneralProc");

  ~G4NeutronGeneralProcess() override;

  G4bool IsApplicable(const G4ParticleDefinition&) override;

  void ProcessDescription(std::ostream& outFile) const override;

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // Store internal tables after initialisation
  G4bool StorePhysicsTable(const G4ParticleDefinition* part,
                           const G4String& directory, G4bool ascii) override;

  // Called before tracking of each new G4Track
  void StartTracking(G4Track*) override;

  // implementation of virtual method, specific for G4NeutronGeneralProcess
  G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition) override;

  // implementation of virtual method, specific for G4NeutronGeneralProcess
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  const G4VProcess* GetCreatorProcess() const override;

  // Temporary method
  const G4String& GetSubProcessName() const;

  // Temporary method
  G4int GetSubProcessSubType() const;

  void SetInelasticProcess(G4HadronicProcess*);
  void SetElasticProcess(G4HadronicProcess*);
  void SetCaptureProcess(G4HadronicProcess*);

  inline const G4VProcess* GetSelectedProcess() const;

  inline void SetTimeLimit(G4double val);

  inline void SetMinEnergyLimit(G4double val);

  // hide copy constructor and assignment operator
  G4NeutronGeneralProcess(G4NeutronGeneralProcess &) = delete;
  G4NeutronGeneralProcess & operator=
  (const G4NeutronGeneralProcess &right) = delete;

protected:

  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                           G4ForceCondition* condition) override;

  inline G4double ComputeGeneralLambda(size_t idxe, size_t idxt);

  inline G4double GetProbability(size_t idxt);

  inline void SelectedProcess(const G4Step& step, G4HadronicProcess* ptr,
                              G4CrossSectionDataStore*);

private:

  // partial cross section
  G4double ComputeCrossSection(G4VCrossSectionDataSet*, const G4Material*,
                               G4double kinEnergy, G4double loge);

  G4VCrossSectionDataSet* InitialisationXS(G4HadronicProcess*);

  // total cross section
  inline void CurrentCrossSection(const G4Track&);

  static G4HadDataHandler* theHandler;
  static const size_t nTables = 5;
  static G4String nameT[nTables];

  G4HadronicProcess* fInelastic = nullptr;
  G4HadronicProcess* fElastic = nullptr;
  G4HadronicProcess* fCapture = nullptr;
  G4HadronicProcess* fSelectedProc = nullptr;

  G4VCrossSectionDataSet* fInelasticXS = nullptr;
  G4VCrossSectionDataSet* fElasticXS = nullptr;
  G4VCrossSectionDataSet* fCaptureXS = nullptr;

  G4CrossSectionDataStore* fXSSInelastic = nullptr;
  G4CrossSectionDataStore* fXSSElastic = nullptr;
  G4CrossSectionDataStore* fXSSCapture = nullptr;
  G4CrossSectionDataStore* fCurrentXSS = nullptr;

  const G4ParticleDefinition* fNeutron;
  const G4Material* fCurrMat = nullptr;

  G4double fMinEnergy;
  G4double fMiddleEnergy;
  G4double fMaxEnergy;
  G4double fTimeLimit;
  G4double fXSFactorInel = 1.0;
  G4double fXSFactorEl = 1.0;
  G4double fCurrE = 0.0;
  G4double fCurrLogE = 0.0;
  G4double fLambda = 0.0;

  // number of bins per decade
  std::size_t nLowE = 100;
  std::size_t nHighE = 10;

  std::size_t idxEnergy = 0;
  std::size_t matIndex = 0;

  G4bool isMaster = true;
  std::vector<G4double> fXsec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double
G4NeutronGeneralProcess::ComputeGeneralLambda(std::size_t idxe, std::size_t idxt)
{
  idxEnergy = idxe;
  return theHandler->GetVector(idxt, matIndex)
    ->LogVectorValue(fCurrE, fCurrLogE);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4NeutronGeneralProcess::GetProbability(std::size_t idxt)
{
  return theHandler->GetVector(idxt, matIndex)
    ->LogVectorValue(fCurrE, fCurrLogE);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void
G4NeutronGeneralProcess::SelectedProcess(const G4Step& step,
                                         G4HadronicProcess* ptr,
                                         G4CrossSectionDataStore* xs)

{
  fSelectedProc = ptr;
  fCurrentXSS = xs;
  step.GetPostStepPoint()->SetProcessDefinedStep(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4VProcess* G4NeutronGeneralProcess::GetSelectedProcess() const
{
  return fSelectedProc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4NeutronGeneralProcess::CurrentCrossSection(const G4Track& track)
{
  G4double energy = track.GetKineticEnergy();
  const G4Material* mat = track.GetMaterial();
  if(mat != fCurrMat || energy != fCurrE) {
    fCurrMat = mat;
    matIndex = mat->GetIndex();
    fCurrE = energy;
    fCurrLogE = track.GetDynamicParticle()->GetLogKineticEnergy();
    fLambda = (energy <= fMiddleEnergy) ? ComputeGeneralLambda(0, 0)
      : ComputeGeneralLambda(1, 3);
    currentInteractionLength = 1.0/fLambda;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4NeutronGeneralProcess::SetTimeLimit(G4double val)
{
  fTimeLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4NeutronGeneralProcess::SetMinEnergyLimit(G4double val)
{
  fMinEnergy = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
