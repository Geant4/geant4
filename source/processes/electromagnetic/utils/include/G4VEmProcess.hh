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
// GEANT4 Class header file
//
//
// File name:     G4VEmProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
//
// Modifications: Vladimir Ivanchenko
//
// Class Description:
//
// It is the base class - EM discrete and rest/discrete process

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
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4EmParameters.hh"
#include "G4EmDataHandler.hh"
#include "G4EmTableType.hh"
#include "G4EmModelManager.hh"
#include "G4EmSecondaryParticleType.hh"

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

  G4VEmProcess(const G4String& name, G4ProcessType type = fElectromagnetic);

  ~G4VEmProcess() override;

  //------------------------------------------------------------------------
  // Virtual methods to be implemented in concrete processes
  //------------------------------------------------------------------------

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override = 0;

  void ProcessDescription(std::ostream& outFile) const override;

protected:

  virtual void StreamProcessInfo(std::ostream&) const {};

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Implementation of virtual methods common to all Discrete processes 
  //------------------------------------------------------------------------

public:

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // Called before tracking of each new G4Track
  void StartTracking(G4Track*) override;

  // implementation of virtual method, specific for G4VEmProcess
  G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) override;

  // implementation of virtual method, specific for G4VEmProcess
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  // Store PhysicsTable in a file.
  // Return false in case of failure at I/O
  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                           G4bool ascii = false) override;

  // Retrieve Physics from a file.
  // (return true if the Physics Table can be build by using file)
  // (return false if the process has no functionality or in case of failure)
  // File name should is constructed as processName+particleName and the
  // should be placed under the directory specified by the argument.
  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                              G4bool ascii) override;

  // allowing check process name
  virtual G4VEmProcess* GetEmProcess(const G4String& name);

  //------------------------------------------------------------------------
  // Specific methods for Discrete EM post step simulation 
  //------------------------------------------------------------------------

  // The main method to access cross section per volume
  inline G4double GetLambda(G4double kinEnergy,
                            const G4MaterialCutsCouple* couple,
                            G4double logKinEnergy);

  // It returns the cross section per volume for energy/material
  G4double GetCrossSection(const G4double kinEnergy,
                           const G4MaterialCutsCouple* couple) override;

  // It returns the cross section of the process per atom
  G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, 
                                      G4double Z, G4double A=0., 
                                      G4double cut=0.0);

  inline G4double MeanFreePath(const G4Track& track);

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
  inline void SetLambdaTable(G4PhysicsTable*);
  inline void SetLambdaTablePrim(G4PhysicsTable*);

  // Integral method type and peak positions
  inline std::vector<G4double>* EnergyOfCrossSectionMax() const;
  inline void SetEnergyOfCrossSectionMax(std::vector<G4double>*);
  inline G4CrossSectionType CrossSectionType() const;
  inline void SetCrossSectionType(G4CrossSectionType val);

  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

  inline const G4ParticleDefinition* Particle() const;
  inline const G4ParticleDefinition* SecondaryParticle() const;

protected:

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models and basic parameters
  //------------------------------------------------------------------------

  // Select model in run time
  inline G4VEmModel* SelectModel(G4double kinEnergy, size_t);

public:

  // Select model by energy and couple index
  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
                                            size_t idxCouple) const;
   
  // Add model for region, smaller value of order defines which
  // model will be selected for a given energy interval  
  void AddEmModel(G4int, G4VEmModel*, const G4Region* region = nullptr);

  // Assign a model to a process local list, to enable the list in run time 
  // the derived process should execute AddEmModel(..) for all such models
  void SetEmModel(G4VEmModel*, G4int index = 0);

  inline G4int NumberOfModels() const;
      
  // return a model from the local list
  inline G4VEmModel* EmModel(size_t index = 0) const;

  // Access to active model
  inline const G4VEmModel* GetCurrentModel() const;

  // Access to models
  inline G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false) const;

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

  inline void SetEmMasterProcess(const G4VEmProcess*);

  inline void SetBuildTableFlag(G4bool val);

  inline void CurrentSetup(const G4MaterialCutsCouple*, G4double energy);

  inline G4bool UseBaseMaterial() const;

  void BuildLambdaTable();

  void StreamInfo(std::ostream& outFile, const G4ParticleDefinition&,
                  G4bool rst=false) const;

  // hide copy constructor and assignment operator
  G4VEmProcess(G4VEmProcess &) = delete;
  G4VEmProcess & operator=(const G4VEmProcess &right) = delete;

  //------------------------------------------------------------------------
  // Other generic methods
  //------------------------------------------------------------------------
  
protected:

  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  inline G4int LambdaBinning() const;

  inline G4double MinKinEnergy() const;

  inline G4double MaxKinEnergy() const;

  // Single scattering parameters
  inline G4double PolarAngleLimit() const;

  inline G4ParticleChangeForGamma* GetParticleChange();

  inline void SetParticle(const G4ParticleDefinition* p);
  
  inline void SetSecondaryParticle(const G4ParticleDefinition* p);

  inline size_t CurrentMaterialCutsCoupleIndex() const;

  inline const G4MaterialCutsCouple* MaterialCutsCouple() const;

  inline G4bool ApplyCuts() const;

  inline G4double GetGammaEnergyCut();

  inline G4double GetElectronEnergyCut();

  inline void SetStartFromNullFlag(G4bool val);

  inline void SetSplineFlag(G4bool val);

  const G4Element* GetTargetElement() const;

  const G4Isotope* GetTargetIsotope() const;

  // these two methods assume that vectors are initilized
  // and idx is within vector length
  inline G4int DensityIndex(G4int idx) const;
  inline G4double DensityFactor(G4int idx) const;

private:

  void PrintWarning(G4String tit, G4double val);

  void ComputeIntegralLambda(G4double kinEnergy, G4double logKinEnergy);

  inline G4double GetLambdaFromTable(G4double kinEnergy);

  inline G4double GetLambdaFromTable(G4double kinEnergy, G4double logKinEnergy);

  inline G4double GetLambdaFromTablePrim(G4double kinEnergy);

  inline G4double GetLambdaFromTablePrim(G4double kinEnergy, G4double logKinEnergy);

  inline G4double GetCurrentLambda(G4double kinEnergy);

  inline G4double GetCurrentLambda(G4double kinEnergy, G4double logKinEnergy);

  inline G4double ComputeCurrentLambda(G4double kinEnergy);

  // ======== pointers =========

  G4EmModelManager*            modelManager = nullptr;
  const G4ParticleDefinition*  particle = nullptr;
  const G4ParticleDefinition*  currentParticle = nullptr;
  const G4ParticleDefinition*  theGamma = nullptr;
  const G4ParticleDefinition*  theElectron = nullptr;
  const G4ParticleDefinition*  thePositron = nullptr;
  const G4ParticleDefinition*  secondaryParticle = nullptr;
  const G4VEmProcess*          masterProc = nullptr;
  G4EmDataHandler*             theData = nullptr;
  G4VEmModel*                  currentModel = nullptr;
  G4LossTableManager*          lManager = nullptr;
  G4EmParameters*              theParameters = nullptr;
  const G4Material*            baseMaterial = nullptr;

  // ======== tables and vectors ========
  G4PhysicsTable*              theLambdaTable = nullptr;
  G4PhysicsTable*              theLambdaTablePrim = nullptr;

  const std::vector<G4double>* theCuts = nullptr;
  const std::vector<G4double>* theCutsGamma = nullptr;
  const std::vector<G4double>* theCutsElectron = nullptr;
  const std::vector<G4double>* theCutsPositron = nullptr;

protected:

  // ======== pointers =========

  const G4MaterialCutsCouple*  currentCouple = nullptr;
  const G4Material*            currentMaterial = nullptr;
  G4EmBiasingManager*          biasManager = nullptr;
  std::vector<G4double>*       theEnergyOfCrossSectionMax = nullptr;

private:

  const std::vector<G4double>* theDensityFactor = nullptr;
  const std::vector<G4int>*    theDensityIdx = nullptr;

  // ======== parameters =========
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double minKinEnergyPrim = DBL_MAX;
  G4double lambdaFactor = 0.8;
  G4double invLambdaFactor;
  G4double biasFactor = 1.0;
  G4double massRatio = 1.0;
  G4double fFactor = 1.0;
  G4double fLambda = 0.0;
  G4double fLambdaEnergy = 0.0;

protected:

  G4double mfpKinEnergy = DBL_MAX;
  G4double preStepKinEnergy = 0.0;
  G4double preStepLogKinEnergy = LOG_EKIN_MIN;
  G4double preStepLambda = 0.0;

private:

  G4CrossSectionType fXSType = fEmNoIntegral;

  G4int numberOfModels = 0;
  G4int nLambdaBins = 84;

protected:

  G4int mainSecondaries = 1;
  G4int secID = _EM;
  G4int fluoID = _Fluorescence; 
  G4int augerID = _AugerElectron;
  G4int biasID = _EM;
  G4int tripletID = _TripletElectron;
  size_t currentCoupleIndex = 0;
  size_t basedCoupleIndex = 0;
  size_t coupleIdxLambda = 0;
  size_t idxLambda = 0;

  G4bool isTheMaster = true;
  G4bool baseMat = false;

private:

  G4bool buildLambdaTable = true;
  G4bool applyCuts = false;
  G4bool startFromNull = false;
  G4bool splineFlag = true;
  G4bool actMinKinEnergy = false;
  G4bool actMaxKinEnergy = false;
  G4bool actBinning = false;
  G4bool isIon = false;
  G4bool biasFlag = false;
  G4bool weightFlag = false;

protected:

  // ======== particle change =========
  std::vector<G4DynamicParticle*> secParticles;
  G4ParticleChangeForGamma fParticleChange;

private:

  // ======== local vectors =========
  std::vector<G4VEmModel*> emModels;

};

// ======== Run time inline methods ================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline size_t G4VEmProcess::CurrentMaterialCutsCoupleIndex() const 
{
  return currentCoupleIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4MaterialCutsCouple* G4VEmProcess::MaterialCutsCouple() const
{
  return currentCouple;
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
    currentCouple = couple;
    baseMaterial = currentMaterial = couple->GetMaterial();
    basedCoupleIndex = currentCoupleIndex = couple->GetIndex();
    fFactor = biasFactor;
    mfpKinEnergy = DBL_MAX;
    if(baseMat) {
      basedCoupleIndex = (*theDensityIdx)[currentCoupleIndex];
      if(nullptr != currentMaterial->GetBaseMaterial())
        baseMaterial = currentMaterial->GetBaseMaterial();
      fFactor *= (*theDensityFactor)[currentCoupleIndex];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VEmProcess::SelectModel(G4double kinEnergy, size_t)
{
  if(1 < numberOfModels) {
    currentModel = modelManager->SelectModel(kinEnergy, currentCoupleIndex);
  }
  currentModel->SetCurrentCouple(currentCouple);
  return currentModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VEmProcess::SelectModelForMaterial(G4double kinEnergy, 
                                                 size_t idxCouple) const
{
  return modelManager->SelectModel(kinEnergy, idxCouple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTable(G4double e)
{
  return ((*theLambdaTable)[basedCoupleIndex])->Value(e, idxLambda);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTable(G4double e, G4double loge)
{
  return ((*theLambdaTable)[basedCoupleIndex])->LogVectorValue(e, loge);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTablePrim(G4double e)
{
  return ((*theLambdaTablePrim)[basedCoupleIndex])->Value(e, idxLambda)/e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTablePrim(G4double e, G4double loge)
{
  return ((*theLambdaTablePrim)[basedCoupleIndex])->LogVectorValue(e, loge)/e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::ComputeCurrentLambda(G4double e)
{
  return currentModel->CrossSectionPerVolume(baseMaterial, currentParticle, e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetCurrentLambda(G4double e)
{
  if(currentCoupleIndex != coupleIdxLambda || fLambdaEnergy != e) {
    coupleIdxLambda = currentCoupleIndex;
    fLambdaEnergy = e;
    if(e >= minKinEnergyPrim) { fLambda = GetLambdaFromTablePrim(e); }
    else if(nullptr != theLambdaTable) { fLambda = GetLambdaFromTable(e); }
    else { fLambda = ComputeCurrentLambda(e); }
    fLambda *= fFactor;
  }
  return fLambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetCurrentLambda(G4double e, G4double loge)
{
  if(currentCoupleIndex != coupleIdxLambda || fLambdaEnergy != e) {
    coupleIdxLambda = currentCoupleIndex;
    fLambdaEnergy = e;
    if(e >= minKinEnergyPrim) { fLambda = GetLambdaFromTablePrim(e, loge); }
    else if(nullptr != theLambdaTable) { fLambda = GetLambdaFromTable(e, loge); }
    else { fLambda = ComputeCurrentLambda(e); }
    fLambda *= fFactor;
  }
  return fLambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEmProcess::CurrentSetup(const G4MaterialCutsCouple* couple, G4double energy)
{
  DefineMaterial(couple);
  SelectModel(energy*massRatio, currentCoupleIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEmProcess::GetLambda(G4double kinEnergy, const G4MaterialCutsCouple* couple,
                        G4double logKinEnergy)
{
  CurrentSetup(couple, kinEnergy);
  return GetCurrentLambda(kinEnergy, logKinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MeanFreePath(const G4Track& track)
{
  const G4double kinEnergy = track.GetKineticEnergy();
  CurrentSetup(track.GetMaterialCutsCouple(), kinEnergy);
  const G4double xs = GetCurrentLambda(kinEnergy,
                             track.GetDynamicParticle()->GetLogKineticEnergy());
  return (0.0 < xs) ? 1.0/xs : DBL_MAX; 
}

// ======== Get/Set inline methods used at initialisation ================

inline G4bool G4VEmProcess::ApplyCuts() const 
{
  return applyCuts;
}

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

inline void G4VEmProcess::SetLambdaTable(G4PhysicsTable* ptr)
{
  theLambdaTable = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetLambdaTablePrim(G4PhysicsTable* ptr)
{
  theLambdaTablePrim = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4double>* G4VEmProcess::EnergyOfCrossSectionMax() const
{
  return theEnergyOfCrossSectionMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEmProcess::SetEnergyOfCrossSectionMax(std::vector<G4double>* ptr)
{
  theEnergyOfCrossSectionMax = ptr;
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

inline void G4VEmProcess::SetCrossSectionType(G4CrossSectionType val)
{
  fXSType = val; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4CrossSectionType G4VEmProcess::CrossSectionType() const
{
  return fXSType;
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmProcess::DensityIndex(G4int idx) const
{
  return (*theDensityIdx)[idx];  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::DensityFactor(G4int idx) const
{
  return (*theDensityFactor)[idx];  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmProcess::UseBaseMaterial() const
{
  return baseMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4VEmModel* G4VEmProcess::GetCurrentModel() const
{
  return currentModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetEmMasterProcess(const G4VEmProcess* ptr)
{ 
  masterProc = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmProcess::NumberOfModels() const
{
  return numberOfModels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmProcess::EmModel(size_t index) const
{
  return (index < emModels.size()) ? emModels[index] : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmProcess::GetModelByIndex(G4int idx, G4bool ver) const
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
