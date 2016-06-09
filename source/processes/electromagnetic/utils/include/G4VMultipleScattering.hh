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
// $Id: G4VMultipleScattering.hh,v 1.29 2005/10/27 11:33:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
// 16-07-03 Update GetRange interface (V.Ivanchenko)
//
//
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles
//
// 26-11-03 bugfix in AlongStepDoIt (L.Urban)
// 25-05-04 add protection against case when range is less than steplimit (V.Ivanchenko)
// 30-06-04 make destructor virtual (V.Ivanchenko)
// 27-08-04 Add InitialiseForRun method (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 15-04-05 optimize internal interfaces (V.Ivanchenko)
// 15-04-05 remove boundary flag (V.Ivanchenko)
// 07-10-05 error in a protection in GetContinuousStepLimit corrected (L.Urban)
// 27-10-05 introduce virtual function MscStepLimitation() (V.Ivanchenko)

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
#include "G4Step.hh"
#include "G4EmModelManager.hh"
#include "G4VEmModel.hh"

class G4ParticleDefinition;
class G4DataVector;
class G4PhysicsTable;
class G4PhysicsVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VMultipleScattering : public G4VContinuousDiscreteProcess
{
public:

  G4VMultipleScattering(const G4String& name = "msc",
                              G4ProcessType type = fElectromagnetic);

  virtual ~G4VMultipleScattering();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;
    // True for all charged particles

  virtual G4double TruePathLengthLimit(const G4Track& track,
                                             G4double& lambda,
                                             G4double currentMinimalStep) = 0;

  virtual void PrintInfo() = 0;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed
  //------------------------------------------------------------------------
public:

  // Initialise for build of tables
  virtual void PreparePhysicsTable(const G4ParticleDefinition&);
  
  // Build physics table during initialisation
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  // set boolean flag steppingAlgorithm
  // ( true/false : standard or 7.1 style process)
  virtual void MscStepLimitation(G4bool algorithm, G4double factor = -1.);

  //------------------------------------------------------------------------
  // Generic methods common to all models
  //------------------------------------------------------------------------

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // The function overloads the corresponding function of the base
  // class.It limits the step near to boundaries only
  // and invokes the method GetContinuousStepLimit at every step.
  G4double AlongStepGetPhysicalInteractionLength(
                                            const G4Track&,
                                                  G4double  previousStepSize,
                                                  G4double  currentMinimalStep,
                                                  G4double& currentSafety,
                                                  G4GPILSelection* selection);

  // Print out of generic class parameters
  void PrintInfoDefinition();

  void SetBinning(G4int nbins);
  G4int Binning() const;
    // Print out of the class parameters

  void SetMinKinEnergy(G4double e);
  G4double MinKinEnergy() const;
    // Print out of the class parameters

  void SetMaxKinEnergy(G4double e);
  G4double MaxKinEnergy() const;

  // Build empty Physics Vector
  G4PhysicsVector* PhysicsVector(const G4MaterialCutsCouple*);

  // Store PhysicsTable in a file.
  // Return false in case of failure at I/O
  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                                 G4bool ascii = false);

  // Retrieve Physics from a file.
  // (return true if the Physics Table can be build by using file)
  // (return false if the process has no functionality or in case of failure)
  // File name should is constructed as processName+particleName and the
  // should be placed under the directory specifed by the argument.
  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                                    G4bool ascii);

  void AddEmModel(G4int, G4VEmModel*, const G4Region* region = 0);

  // This method does not used for tracking, it is intended only for tests
  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimalStep,
                                     G4double& currentSafety);

  G4bool LateralDisplasmentFlag() const;
  void SetLateralDisplasmentFlag(G4bool val);
     // lateral displacement to be/not to be computed

  void SetBuildLambdaTable(G4bool val);

  const G4PhysicsTable* LambdaTable() const;

  G4VEmModel* SelectModelForMaterial(G4double kinEnergy, size_t& idxRegion) const;

  // Define particle definition
  const G4ParticleDefinition* Particle() const;
  void SetParticle(const G4ParticleDefinition*);

protected:

  // This method is used for tracking, it returns mean free path value
  G4double GetMeanFreePath(const G4Track& track,
                                 G4double,
                                 G4ForceCondition* condition);

  G4double GetLambda(const G4ParticleDefinition* p, G4double& kineticEnergy);

  // This method is used for tracking, it returns step limit
  G4double GetContinuousStepLimit(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimalStep,
                                        G4double& currentSafety);

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

  G4ParticleChangeForMSC      fParticleChange;
  G4EmModelManager*           modelManager;
  G4VEmModel*                 currentModel;

  // tables and vectors
  G4PhysicsTable*             theLambdaTable;

  // cash
  const G4ParticleDefinition* firstParticle;
  const G4ParticleDefinition* currentParticle;
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
  G4bool                      latDisplasment;
  G4bool                      buildLambdaTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterialIndex = couple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::GetMeanFreePath(const G4Track&,
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
  currentRange = G4LossTableManager::Instance()->GetTrancatedRange(p,e,currentCouple);
  // the next line was in error 
  // if(currentRange < currentMinimalStep) currentRange = currentMinimalStep;
  //  the condition/protection correctly should be
  if(currentRange < currentMinimalStep) currentMinimalStep = currentRange;
  truePathLength = TruePathLengthLimit(track,lambda0,currentMinimalStep);
  //G4cout << "StepLimit: tpl= " << truePathLength << " lambda0= "
  //       << lambda0 << " range= " << currentRange
  //       << " currentMinStep= " << currentMinimalStep << G4endl;
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
    x = currentModel->CrossSection(currentCouple,p,e);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4VMultipleScattering::AlongStepDoIt(
                                                        const G4Track&,
                                                        const G4Step& step)
{
  G4double geomStepLength = step.GetStepLength();
  if((geomStepLength == geomPathLength) && (truePathLength <= currentRange))
     trueStepLength = truePathLength;
  else
     trueStepLength = currentModel->TrueStepLength(geomStepLength);
  fParticleChange.ProposeTrueStepLength(trueStepLength);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4VMultipleScattering::PostStepDoIt(const G4Track& track,
							      const G4Step& step)
{
  fParticleChange.Initialize(track);
  currentModel->SampleSecondaries(currentCouple,track.GetDynamicParticle(),
		    step.GetStepLength(),step.GetPostStepPoint()->GetSafety());
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SelectModel(G4double& kinEnergy)
{
  currentModel = modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VMultipleScattering::SelectModelForMaterial(
                                           G4double kinEnergy, size_t& idxRegion) const
{
  return modelManager->SelectModel(kinEnergy, idxRegion);
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

inline  void G4VMultipleScattering::SetBuildLambdaTable(G4bool val)
{
  buildLambdaTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  const G4ParticleDefinition* G4VMultipleScattering::Particle() const
{
  return currentParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4PhysicsTable* G4VMultipleScattering::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
