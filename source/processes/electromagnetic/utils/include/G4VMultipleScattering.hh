//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VMultipleScattering
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 12.03.2002
//
// Modifications:
//
//
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles

// -------------------------------------------------------------------
//

#ifndef G4VMultipleScattering_h
#define G4VMultipleScattering_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4LossTableManager.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4DataVector;
class G4Navigator;
class G4PhysicsTable;
class G4PhysicsVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VMultipleScattering : public G4VContinuousDiscreteProcess
{
public:

  G4VMultipleScattering(const G4String& name = "msc",
                              G4ProcessType type = fElectromagnetic);

 ~G4VMultipleScattering();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;
    // True for all charged particles

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    // Build physics table during initialisation

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double AlongStepGetPhysicalInteractionLength(
                                            const G4Track&,
                                                  G4double  previousStepSize,
                                                  G4double  currentMinimalStep,
                                                  G4double& currentSafety,
                                                  G4GPILSelection* selection);
     // The function overloads the corresponding function of the base
     // class.It limits the step near to boundaries only
     // and invokes the method GetContinuousStepLimit at every step.

  virtual void PrintInfoDefinition() const;
  // Print out of the class parameters

  void SetBinning(G4int nbins);
  G4int Binning() const;
    // Print out of the class parameters

  void SetMinKinEnergy(G4double e);
  G4double MinKinEnergy() const;
    // Print out of the class parameters

  void SetMaxKinEnergy(G4double e);
  G4double MaxKinEnergy() const;
    // Print out of the class parameters

  G4bool StorePhysicsTable(G4ParticleDefinition*,
                     const G4String& directory,
                           G4bool ascii = false);
    // Store PhysicsTable in a file.
    // Return false in case of failure at I/O

  G4bool RetrievePhysicsTable(G4ParticleDefinition*,
                        const G4String& directory,
                              G4bool ascii);
    // Retrieve Physics from a file.
    // (return true if the Physics Table can be build by using file)
    // (return false if the process has no functionality or in case of failure)
    // File name should is constructed as processName+particleName and the
    // should be placed under the directory specifed by the argument.

  void AddEmModel(G4int, G4VEmModel*, const G4Region* region = 0);

  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimalStep,
                                     G4double& currentSafety);
  // This method does not used for tracking, it is intended only for tests

  G4bool LateralDisplasmentFlag() const;
  void SetLateralDisplasmentFlag(G4bool val);
     // lateral displacement to be/not to be computed

  G4bool BoundaryAlgorithmFlag() const;
  void SetBoundary(G4bool val);
     // boundary algorith is/isnt active

  void SetBuildLambdaTable(G4bool val);

  virtual G4double TruePathLengthLimit(const G4Track& track,
                                             G4double& lambda,
                                             G4double currentMinimalStep) = 0;


protected:

  virtual void InitialiseProcess(const G4ParticleDefinition&) = 0;

  G4double GetMeanFreePath(const G4Track& track,
                                 G4double,
                                 G4ForceCondition* condition);
  // This method is used for tracking, it returns mean free path value

  G4double GetLambda(const G4ParticleDefinition* p, G4double& kineticEnergy);

  G4double GetContinuousStepLimit(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimalStep,
                                        G4double& currentSafety);
  // This method is used for tracking, it returns step limit

  virtual G4PhysicsVector* PhysicsVector(const G4MaterialCutsCouple*);
  // Build empty Physics Vector

  void SelectModel(G4double& kinEnergy);
  // Select concrete model

  size_t CurrentMaterialCutsCoupleIndex() const {return currentMaterialIndex;};
  // Return current index

  G4double CurrentRange() const {return currentRange;};

private:

  void DefineMaterial(const G4MaterialCutsCouple* couple);
  // define current material

  // hide  assignment operator

  G4VMultipleScattering(G4VMultipleScattering &);
  G4VMultipleScattering & operator=(const G4VMultipleScattering &right);

// =====================================================================

private:

  G4ParticleChangeForMSC      fParticleChange;
  G4EmModelManager*           modelManager;
  G4Navigator*                navigator;
  G4VEmModel*                 currentModel;

  // tables and vectors
  G4PhysicsTable*             theLambdaTable;

  // cash
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;

  G4int                       nBins;

  G4double                    minKinEnergy;
  G4double                    maxKinEnergy;

  G4double                    trueStepLength;
  G4double                    truePathLength;
  G4double                    geomPathLength;
  G4double                    lambda0;
  G4double                    currentRange;

  G4GPILSelection             valueGPILSelectionMSC;
  G4bool                      boundary;
  G4bool                      latDisplasment;
  G4bool                      buildLambdaTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::GetMeanFreePath(const G4Track& track,
                                                             G4double,
                                                             G4ForceCondition* cond)
{
  *cond = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMultipleScattering::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4double currentMinimalStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection)
{
  // get Step limit proposed by the process
  valueGPILSelectionMSC = NotCandidateForSelection;
  G4double steplength = GetContinuousStepLimit(track,previousStepSize,
                                              currentMinimalStep,currentSafety);
  // set return value for G4GPILSelection
  *selection = valueGPILSelectionMSC;
  return  steplength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::GetContinuousStepLimit(
                                          const G4Track& track,
                                                G4double,
                                                G4double currentMinimalStep,
                                                G4double&)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  G4double e = track.GetKineticEnergy();
  SelectModel(e);
  const G4ParticleDefinition* p = track.GetDefinition();
  lambda0 = GetLambda(p, e);
  currentRange = G4LossTableManager::Instance()->GetRange(p,e,currentCouple);
  truePathLength = TruePathLengthLimit(track,lambda0,currentMinimalStep);
  //G4cout << "StepLimit: tpl= " << truePathLength << " lambda0= "
  //       << lambda0 << " range= " << currentRange
  //	 << " currentMinStep= " << currentMinimalStep << G4endl;
  if (truePathLength < currentMinimalStep) valueGPILSelectionMSC = CandidateForSelection;
  geomPathLength = currentModel->GeomPathLength(theLambdaTable,currentCouple,
           p,e,lambda0,currentRange,truePathLength);
  if(geomPathLength > lambda0) geomPathLength = lambda0;
  return geomPathLength;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::ContinuousStepLimit(
                                 const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  return GetContinuousStepLimit(track,previousStepSize,currentMinimalStep,
                                      currentSafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::GetLambda(const G4ParticleDefinition* p, G4double& e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);

  } else {
    x = currentModel->CrossSection(currentMaterial,p,e,0.0,1.0);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4VMultipleScattering::AlongStepDoIt(
                                                        const G4Track& track,
                                                        const G4Step& step)
{
  fParticleChange.Initialize(track);
  G4double geomStepLength = step.GetStepLength();
  if (geomStepLength == geomPathLength) trueStepLength = truePathLength;
  else    trueStepLength = currentModel->TrueStepLength(geomStepLength);
  fParticleChange.SetTrueStepLength(trueStepLength);
//  G4cout << "AlongStep: trueLength= " << trueStepLength << " geomLength= "<< geomStepLength << " zlast= " << geomPathLength << G4endl;
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SelectModel(G4double& kinEnergy)
{
  currentModel = modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetBinning(G4int nbins)
{
  nBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VMultipleScattering::Binning() const
{
  return nBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4bool G4VMultipleScattering::LateralDisplasmentFlag() const
{
  return latDisplasment;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetLateralDisplasmentFlag(G4bool val)
{
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4bool G4VMultipleScattering::BoundaryAlgorithmFlag() const
{
  return boundary;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetBoundary(G4bool val)
{
  boundary = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetBuildLambdaTable(G4bool val)
{
  buildLambdaTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
