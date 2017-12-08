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
// $Id: G4VEmProcess.hh 106983 2017-10-31 09:08:30Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VEmProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
//
// Modifications:
// 30-06-04 make destructor virtual (V.Ivanchenko)
// 09-08-04 optimise integral option (V.Ivanchenko)
// 11-08-04 add protected methods to access cuts (V.Ivanchenko)
// 09-09-04 Bug fix for the integral mode with 2 peaks (V.Ivanchneko)
// 16-09-04 Add flag for LambdaTable and method RecalculateLambda (VI)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
// 09-05-05 Fix problem in logic when path boundary between materials (VI)
// 11-01-06 add A to parameters of ComputeCrossSectionPerAtom (VI)
// 01-02-06 put default value A=0. to keep compatibility with v5.2 (mma)
// 13-05-06 Add method to access model by index (V.Ivanchenko)
// 12-09-06 add SetModel() (mma)
// 25-09-07 More accurate handling zero xsect in 
//          PostStepGetPhysicalInteractionLength (V.Ivanchenko)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
// 15-07-08 Reorder class members for further multi-thread development (VI)
// 17-02-10 Added pointer currentParticle (VI)
//
// Class Description:
//
// It is the unified Discrete process

// -------------------------------------------------------------------
//

#ifndef G4VEmProcess_h
#define G4VEmProcess_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4EmParameters.hh"

class G4Step;
class G4VEmModel;
class G4DataVector;
class G4VParticleChange;
class G4PhysicsTable;
class G4PhysicsVector;
class G4EmBiasingManager;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEmProcess : public G4VDiscreteProcess
{
public:

  G4VEmProcess(const G4String& name,
               G4ProcessType type = fElectromagnetic);

  virtual ~G4VEmProcess();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented in concrete processes
  //------------------------------------------------------------------------

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override = 0;

  // obsolete
  virtual void PrintInfo() {};

  virtual void ProcessDescription(std::ostream& outFile) const override;

protected:

  virtual void StreamProcessInfo(std::ostream&, G4String) const {};

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Method with standard implementation; may be overwritten if needed
  //------------------------------------------------------------------------

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*);

  //------------------------------------------------------------------------
  // Implementation of virtual methods common to all Discrete processes 
  //------------------------------------------------------------------------

public:

  // Initialise for build of tables
  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // Build physics table during initialisation
  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // Called before tracking of each new G4Track
  virtual void StartTracking(G4Track*) override;

  // implementation of virtual method, specific for G4VEmProcess
  virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) override;

  // implementation of virtual method, specific for G4VEmProcess
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, 
                                          const G4Step&) override;

  // Store PhysicsTable in a file.
  // Return false in case of failure at I/O
  virtual G4bool StorePhysicsTable(const G4ParticleDefinition*,
                                   const G4String& directory,
                                   G4bool ascii = false) override;

  // Retrieve Physics from a file.
  // (return true if the Physics Table can be build by using file)
  // (return false if the process has no functionality or in case of failure)
  // File name should is constructed as processName+particleName and the
  // should be placed under the directory specifed by the argument.
  virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                                      const G4String& directory,
                                      G4bool ascii) override;

  //------------------------------------------------------------------------
  // Specific methods for Discrete EM post step simulation 
  //------------------------------------------------------------------------

  // It returns the cross section per volume for energy/ material
  G4double CrossSectionPerVolume(G4double kineticEnergy,
                                 const G4MaterialCutsCouple* couple);

  // It returns the cross section of the process per atom
  G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, 
                                      G4double Z, G4double A=0., 
                                      G4double cut=0.0);

  G4double MeanFreePath(const G4Track& track);

  // It returns cross section per volume
  inline G4double GetLambda(G4double& kinEnergy, 
                            const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Specific methods to build and access Physics Tables
  //------------------------------------------------------------------------

  // Binning for lambda table
  void SetLambdaBinning(G4int nbins);

  // Min kinetic energy for tables
  void SetMinKinEnergy(G4double e);

  // Min kinetic energy for high energy table
  void SetMinKinEnergyPrim(G4double e);

  // Max kinetic energy for tables
  void SetMaxKinEnergy(G4double e);

  // Cross section table pointers
  inline G4PhysicsTable* LambdaTable() const;
  inline G4PhysicsTable* LambdaTablePrim() const;

  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

  inline const G4ParticleDefinition* Particle() const;
  inline const G4ParticleDefinition* SecondaryParticle() const;

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models and basic parameters
  //------------------------------------------------------------------------

protected:
  // Select model in run time
  inline G4VEmModel* SelectModel(G4double& kinEnergy, size_t index); 

public:
  // Select model by energy and region index
  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
                                            size_t& idxRegion) const;
   
  // Add model for region, smaller value of order defines which
  // model will be selected for a given energy interval  
  void AddEmModel(G4int, G4VEmModel*, const G4Region* region = nullptr);

  // Assign a model to a process local list, to enable the list in run time 
  // the derived process should execute AddEmModel(..) for all such models
  void SetEmModel(G4VEmModel*, G4int index = 0);
      
  // return a model from the local list
  G4VEmModel* EmModel(size_t index = 0) const;

  // Define new energy range for the model identified by the name
  void UpdateEmModel(const G4String&, G4double, G4double);

  // Access to models
  G4int GetNumberOfModels() const;
  G4int GetNumberOfRegionModels(size_t couple_index) const;
  G4VEmModel* GetRegionModel(G4int idx, size_t couple_index) const;
  G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false) const;

  // Access to active model
  inline const G4VEmModel* GetCurrentModel() const;

  // Access to the current G4Element
  const G4Element* GetCurrentElement() const;

  // Biasing parameters
  void SetCrossSectionBiasingFactor(G4double f, G4bool flag = true);
  inline G4double CrossSectionBiasingFactor() const;

  // Activate forced interaction
  void ActivateForcedInteraction(G4double length = 0.0, 
                                 const G4String& r = "",
                                 G4bool flag = true);

  void ActivateSecondaryBiasing(const G4String& region, G4double factor,
          G4double energyLimit);
          
  inline void SetIntegral(G4bool val);

  inline void SetBuildTableFlag(G4bool val);

  //------------------------------------------------------------------------
  // Other generic methods
  //------------------------------------------------------------------------
  
protected:

  virtual G4double GetMeanFreePath(const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  inline G4int LambdaBinning() const;

  inline G4double MinKinEnergy() const;

  inline G4double MaxKinEnergy() const;

  // Single scattering parameters
  inline G4double PolarAngleLimit() const;

  inline G4bool IsIntegral() const;

  inline G4double RecalculateLambda(G4double kinEnergy,
                                     const G4MaterialCutsCouple* couple);

  inline G4ParticleChangeForGamma* GetParticleChange();

  inline void SetParticle(const G4ParticleDefinition* p);
  
  inline void SetSecondaryParticle(const G4ParticleDefinition* p);

  inline size_t CurrentMaterialCutsCoupleIndex() const;

  inline G4double GetGammaEnergyCut();

  inline G4double GetElectronEnergyCut();

  inline void SetStartFromNullFlag(G4bool val);

  inline void SetSplineFlag(G4bool val);

  inline const G4Element* GetTargetElement() const;

  inline const G4Isotope* GetTargetIsotope() const;

private:

  void Clear();

  void BuildLambdaTable();

  void StreamInfo(std::ostream& outFile, const G4ParticleDefinition&,
                  G4String endOfLine=G4String("\n")) const;

  void FindLambdaMax();

  void PrintWarning(G4String tit, G4double val);

  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  inline void ComputeIntegralLambda(G4double kinEnergy);

  inline G4double GetLambdaFromTable(G4double kinEnergy);

  inline G4double GetLambdaFromTablePrim(G4double kinEnergy);

  inline G4double GetCurrentLambda(G4double kinEnergy);

  inline G4double ComputeCurrentLambda(G4double kinEnergy);

  // hide copy constructor and assignment operator
  G4VEmProcess(G4VEmProcess &) = delete;
  G4VEmProcess & operator=(const G4VEmProcess &right) = delete;

  // ======== Parameters of the class fixed at construction =========

  G4LossTableManager*          lManager;
  G4EmParameters*              theParameters;  
  G4EmModelManager*            modelManager;
  G4EmBiasingManager*          biasManager;
  const G4ParticleDefinition*  theGamma;
  const G4ParticleDefinition*  theElectron;
  const G4ParticleDefinition*  thePositron;
  const G4ParticleDefinition*  secondaryParticle;

  G4bool                       buildLambdaTable;

  // ======== Parameters of the class fixed at initialisation =======

  std::vector<G4VEmModel*>     emModels;
  G4int                        numberOfModels;

  // tables and vectors
  G4PhysicsTable*              theLambdaTable;
  G4PhysicsTable*              theLambdaTablePrim;
  std::vector<G4double>        theEnergyOfCrossSectionMax;
  std::vector<G4double>        theCrossSectionMax;

  size_t                       idxLambda;
  size_t                       idxLambdaPrim;

  const std::vector<G4double>* theCuts;
  const std::vector<G4double>* theCutsGamma;
  const std::vector<G4double>* theCutsElectron;
  const std::vector<G4double>* theCutsPositron;
  const std::vector<G4double>* theDensityFactor;
  const std::vector<G4int>*    theDensityIdx;

  G4int                        nLambdaBins;

  G4double                     minKinEnergy;
  G4double                     minKinEnergyPrim;
  G4double                     maxKinEnergy;
  G4double                     lambdaFactor;
  G4double                     biasFactor;
  G4double                     massRatio;

  G4bool                       integral;
  G4bool                       applyCuts;
  G4bool                       startFromNull;
  G4bool                       splineFlag;
  G4bool                       actMinKinEnergy;
  G4bool                       actMaxKinEnergy;
  G4bool                       actBinning;
  G4bool                       actSpline;
  G4bool                       isIon;

  // ======== Cashed values - may be state dependent ================

protected:

  G4ParticleChangeForGamma     fParticleChange;

private:

  std::vector<G4DynamicParticle*> secParticles;

  G4VEmModel*                  currentModel;  

  const G4ParticleDefinition*  particle;
  const G4ParticleDefinition*  currentParticle;

  // cache
  const G4Material*            baseMaterial;
  const G4Material*            currentMaterial;
  const G4MaterialCutsCouple*  currentCouple;
  size_t                       currentCoupleIndex;
  size_t                       basedCoupleIndex;

  G4double                     mfpKinEnergy;
  G4double                     preStepKinEnergy;
  G4double                     preStepLambda;
  G4double                     fFactor;
  G4bool                       biasFlag;
  G4bool                       weightFlag;

  G4int                        mainSecondaries;
  G4int                        secID;  
  G4int                        fluoID;  
  G4int                        augerID;  
  G4int                        biasID;  
};

// ======== Run time inline methods ================

inline size_t G4VEmProcess::CurrentMaterialCutsCoupleIndex() const 
{
  return currentCoupleIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetGammaEnergyCut()
{
  return (*theCutsGamma)[currentCoupleIndex];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetElectronEnergyCut()
{
  return (*theCutsElectron)[currentCoupleIndex];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    baseMaterial = currentMaterial->GetBaseMaterial();
    currentCoupleIndex = couple->GetIndex();
    basedCoupleIndex   = (*theDensityIdx)[currentCoupleIndex];
    fFactor = biasFactor*(*theDensityFactor)[currentCoupleIndex];
    if(!baseMaterial) { baseMaterial = currentMaterial; }
    mfpKinEnergy = DBL_MAX;
    idxLambda = idxLambdaPrim = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VEmProcess::SelectModel(G4double& kinEnergy, size_t index)
{
  if(1 < numberOfModels) {
    currentModel = modelManager->SelectModel(kinEnergy, index);
  }
  currentModel->SetCurrentCouple(currentCouple);
  return currentModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VEmProcess::SelectModelForMaterial(G4double kinEnergy, 
                                                 size_t& idxRegion) const
{
  return modelManager->SelectModel(kinEnergy, idxRegion);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTable(G4double e)
{
  return ((*theLambdaTable)[basedCoupleIndex])->Value(e, idxLambda);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTablePrim(G4double e)
{
  return ((*theLambdaTablePrim)[basedCoupleIndex])->Value(e, idxLambdaPrim)/e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::ComputeCurrentLambda(G4double e)
{
  return currentModel->CrossSectionPerVolume(
         baseMaterial,currentParticle, e,(*theCuts)[currentCoupleIndex]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetCurrentLambda(G4double e)
{
  G4double x;
  if(e >= minKinEnergyPrim) { x = GetLambdaFromTablePrim(e); }
  else if(theLambdaTable)   { x = GetLambdaFromTable(e); }
  else                      { x = ComputeCurrentLambda(e); }
  return fFactor*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEmProcess::GetLambda(G4double& kinEnergy, 
                        const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  SelectModel(kinEnergy, currentCoupleIndex);
  return GetCurrentLambda(kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEmProcess::RecalculateLambda(G4double e, const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  SelectModel(e, currentCoupleIndex);
  return fFactor*ComputeCurrentLambda(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::ComputeIntegralLambda(G4double e)
{
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentCoupleIndex];
  if (e <= mfpKinEnergy) {
    preStepLambda = GetCurrentLambda(e);

  } else {
    G4double e1 = e*lambdaFactor;
    if(e1 > mfpKinEnergy) {
      preStepLambda = GetCurrentLambda(e);
      G4double preStepLambda1 = GetCurrentLambda(e1);
      if(preStepLambda1 > preStepLambda) {
        mfpKinEnergy = e1;
        preStepLambda = preStepLambda1;
      }
    } else {
      preStepLambda = fFactor*theCrossSectionMax[currentCoupleIndex];
    }
  }
}

// ======== Get/Set inline methods used at initialisation ================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmProcess::LambdaBinning() const
{
  return nLambdaBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::PolarAngleLimit() const
{
  return theParameters->MscThetaLimit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::CrossSectionBiasingFactor() const
{
  return biasFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEmProcess::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEmProcess::LambdaTablePrim() const
{
  return theLambdaTablePrim;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEmProcess::Particle() const
{
  return particle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEmProcess::SecondaryParticle() const
{
  return secondaryParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetIntegral(G4bool val)
{
  integral = val; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmProcess::IsIntegral() const
{
  return integral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetBuildTableFlag(G4bool val)
{
  buildLambdaTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4ParticleChangeForGamma* G4VEmProcess::GetParticleChange()
{
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  currentParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetStartFromNullFlag(G4bool val)
{
  startFromNull = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetSplineFlag(G4bool val)
{
  splineFlag = val;
  actSpline = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Element* G4VEmProcess::GetTargetElement() const
{
  return currentModel->GetCurrentElement();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Isotope* G4VEmProcess::GetTargetIsotope() const
{
  return currentModel->GetCurrentIsotope();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4VEmModel* G4VEmProcess::GetCurrentModel() const
{
  return currentModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
