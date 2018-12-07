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
// File name:     G4GammaGeneralProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 19.07.2018
//
// Modifications:
//
// Class Description:
//
// It is the gamma super process

// -------------------------------------------------------------------
//

#ifndef G4GammaGeneralProcess_h
#define G4GammaGeneralProcess_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmProcess.hh"
#include "globals.hh"
#include "G4EmDataHandler.hh"

class G4Step;
class G4Track;
class G4ParticleDefinition;
class G4VParticleChange;
class G4GammaConversionToMuons;
class G4HadronicProcess;
class G4LossTableManager;
class G4MaterialCutsCouple;
class G4EmParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4GammaGeneralProcess : public G4VEmProcess
{
public:

  explicit G4GammaGeneralProcess();

  virtual ~G4GammaGeneralProcess();

  G4bool IsApplicable(const G4ParticleDefinition&) final;

  void AddEmProcess(G4VEmProcess*);

  void AddMMProcess(G4GammaConversionToMuons*);

  void AddHadProcess(G4HadronicProcess*);

  void ProcessDescription(std::ostream& outFile) const final;

protected:

  void InitialiseProcess(const G4ParticleDefinition*) final;

public:

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&) final;

  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  // Called before tracking of each new G4Track
  void StartTracking(G4Track*) final;
  
  // implementation of virtual method, specific for G4GammaGeneralProcess
  G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) final;

  // implementation of virtual method, specific for G4GammaGeneralProcess
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) final;
  
  // Store PhysicsTable in a file.
  // Return false in case of failure at I/O
  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                           G4bool ascii = false) final;

  // Retrieve Physics from a file.
  // (return true if the Physics Table can be build by using file)
  // (return false if the process has no functionality or in case of failure)
  // File name should is constructed as processName+particleName and the
  // should be placed under the directory specifed by the argument.
  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                              G4bool ascii) final;

  // It returns the cross section per volume for energy/ material
  G4double TotalCrossSectionPerVolume(G4double energy,
				      const G4MaterialCutsCouple* couple);

  const G4String& GetProcessName() const;

  G4int GetProcessSubType() const;

protected:

  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                           G4ForceCondition* condition) final;

private:

  inline G4double ComputeGeneralLambda(size_t idxe, size_t idxt, 
                                     size_t& idx, G4double e);

  inline G4double GetProbability(size_t idxt, size_t& idx);

  inline void SelectedProcess(const G4Track& track, G4VProcess* ptr);

  inline G4VParticleChange* SampleSecondaries(const G4Track&, const G4Step&,
					      G4VEmProcess*);

  G4VParticleChange* SampleSecondaries(const G4Track&, const G4Step&,
				       G4HadronicProcess*);

  // hide copy constructor and assignment operator
  G4GammaGeneralProcess(G4GammaGeneralProcess &) = delete;
  G4GammaGeneralProcess & operator=(const G4GammaGeneralProcess &right) = delete;

  static G4EmDataHandler*      theHandler;
  static const size_t          nTables = 15;
  static G4bool                theT[nTables];
  static G4String              nameT[nTables];

  G4VEmProcess*                thePhotoElectric;
  G4VEmProcess*                theCompton;
  G4VEmProcess*                theConversionEE;
  G4VEmProcess*                theRayleigh;
  G4HadronicProcess*           theGammaNuclear;
  G4GammaConversionToMuons*    theConversionMM;
  G4VProcess*                  selectedProc;

  G4double                     minPEEnergy;
  G4double                     minEEEnergy;
  G4double                     minMMEnergy;
  G4double                     peLambda;

  size_t                       nLowE;
  size_t                       nHighE;
  size_t                       idxEnergy;
  size_t                       idx0;
  size_t                       idx1;
  size_t                       idx2;
  size_t                       idx3;
  G4bool                       splineFlag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
inline G4double
G4GammaGeneralProcess::ComputeGeneralLambda(size_t idxe, size_t idxt, 
				            size_t& idx, G4double e)
{
  idxEnergy = idxe;
  return theHandler->GetVector(idxt, currentCoupleIndex)->Value(e, idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4GammaGeneralProcess::GetProbability(size_t idxt, size_t& idx)
{
  return (theT[idxt]) ? theHandler->GetVector(idxt, 
          currentCoupleIndex)->Value(preStepKinEnergy, idx) : 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4GammaGeneralProcess::SelectedProcess(const G4Track& track, G4VProcess* ptr)
{
  selectedProc = ptr;
  const_cast<G4Track*>(&track)->SetCreatorProcess(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* 
G4GammaGeneralProcess::SampleSecondaries(const G4Track& track, const G4Step& step,
	 			         G4VEmProcess* proc)
{
  proc->CurrentSetup(currentCouple,preStepKinEnergy);
  SelectedProcess(track, proc);
  return proc->PostStepDoIt(track, step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
