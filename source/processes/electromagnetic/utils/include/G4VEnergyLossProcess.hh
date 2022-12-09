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
// File name:     G4VEnergyLossProcess
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 03.01.2002
//
// Modifications: Vladimir Ivanchenko
//
// Class Description:
//
// It is the unified energy loss process it calculates the continuous
// energy loss for charged particles using a set of Energy Loss
// models valid for different energy regions. There are a possibility
// to create and access to dE/dx and range tables, or to calculate
// that information on fly.

// -------------------------------------------------------------------
//

#ifndef G4VEnergyLossProcess_h
#define G4VEnergyLossProcess_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4EmTableType.hh"
#include "G4EmSecondaryParticleType.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"

class G4Step;
class G4ParticleDefinition;
class G4EmParameters;
class G4VEmModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4Region;
class G4SafetyHelper;
class G4VAtomDeexcitation;
class G4VSubCutProducer;
class G4EmBiasingManager;
class G4LossTableManager;
class G4EmDataHandler;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEnergyLossProcess : public G4VContinuousDiscreteProcess
{
public:

  G4VEnergyLossProcess(const G4String& name = "EnergyLoss",
                       G4ProcessType type = fElectromagnetic);

  ~G4VEnergyLossProcess() override;

  //------------------------------------------------------------------------
  // Virtual methods to be implemented in concrete processes
  //------------------------------------------------------------------------

protected:

  // description of specific process parameters
  virtual void StreamProcessInfo(std::ostream&) const {};

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*) = 0;

public:

  // used as low energy limit LambdaTable
  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut);

  // print documentation in html format
  void ProcessDescription(std::ostream& outFile) const override;

  // prepare all tables
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // build all tables
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // build a table
  G4PhysicsTable* BuildDEDXTable(G4EmTableType tType = fRestricted);

  // build a table
  G4PhysicsTable* BuildLambdaTable(G4EmTableType tType = fRestricted);

  // Called before tracking of each new G4Track
  void StartTracking(G4Track*) override;

  // Step limit from AlongStep 
  G4double AlongStepGetPhysicalInteractionLength(
                                  const G4Track&,
                                  G4double  previousStepSize,
                                  G4double  currentMinimumStep,
                                  G4double& currentSafety,
                                  G4GPILSelection* selection) override;

  // Step limit from cross section
  G4double PostStepGetPhysicalInteractionLength(
                                  const G4Track& track,
                                  G4double previousStepSize,
                                  G4ForceCondition* condition) override;

  // AlongStep computations
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  // PostStep sampling of secondaries
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  // Store all PhysicsTable in files.
  // Return false in case of any fatal failure at I/O  
  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                           G4bool ascii = false) override;

  // Retrieve all Physics from a files.
  // Return true if all the Physics Table are built.
  // Return false if any fatal failure. 
  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                              G4bool ascii) override;

private:

  // summary printout after initialisation
  void StreamInfo(std::ostream& out, const G4ParticleDefinition& part,
                  G4bool rst=false) const;

  //------------------------------------------------------------------------
  // Public interface to cross section, mfp and sampling of fluctuations
  // These methods are not used in run time
  //------------------------------------------------------------------------

public:

  // access to dispersion of restricted energy loss
  G4double GetDEDXDispersion(const G4MaterialCutsCouple *couple,
                             const G4DynamicParticle* dp,
                             G4double length);

  // Access to cross section table
  G4double CrossSectionPerVolume(G4double kineticEnergy,
                                 const G4MaterialCutsCouple* couple);
  G4double CrossSectionPerVolume(G4double kineticEnergy,
                                 const G4MaterialCutsCouple* couple,
                                 G4double logKineticEnergy);

  // access to cross section
  G4double MeanFreePath(const G4Track& track);

  // access to step limit
  G4double ContinuousStepLimit(const G4Track& track,
                               G4double previousStepSize,
                               G4double currentMinimumStep,
                               G4double& currentSafety);

protected:

  // implementation of the pure virtual method
  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  // implementation of the pure virtual method
  G4double GetContinuousStepLimit(const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimumStep,
                                  G4double& currentSafety) override;

  // creation of an empty vector for cross sections for derived processes
  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*, 
                                       G4double cut);

  inline std::size_t CurrentMaterialCutsCoupleIndex() const;

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  // Select model in run time
  inline void SelectModel(G4double kinEnergy);

public:
  // Select model by energy and couple index
  // Not for run time processing
  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
                                            std::size_t& idxCouple) const;

  // Add EM model coupled with fluctuation model for region, smaller value 
  // of order defines which pair of models will be selected for a given 
  // energy interval  
  void AddEmModel(G4int, G4VEmModel*, 
                  G4VEmFluctuationModel* fluc = nullptr,
                  const G4Region* region = nullptr);

  // Assign a model to a process local list, to enable the list in run time 
  // the derived process should execute AddEmModel(..) for all such models
  void SetEmModel(G4VEmModel*, G4int index=0);

  // Access to models
  inline std::size_t NumberOfModels() const;
  
  // Return a model from the local list
  inline G4VEmModel* EmModel(std::size_t index=0) const;
  
  // Access to models from G4EmModelManager list
  inline G4VEmModel* GetModelByIndex(std::size_t idx = 0, G4bool ver = false) const;

  // Assign a fluctuation model to a process
  inline void SetFluctModel(G4VEmFluctuationModel*);
  
  // Return the assigned fluctuation model
  inline G4VEmFluctuationModel* FluctModel() const;
    
  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

protected:
  inline void SetParticle(const G4ParticleDefinition* p);
  inline void SetSecondaryParticle(const G4ParticleDefinition* p);

public:
  inline void SetBaseParticle(const G4ParticleDefinition* p);
  inline const G4ParticleDefinition* Particle() const;
  inline const G4ParticleDefinition* BaseParticle() const;
  inline const G4ParticleDefinition* SecondaryParticle() const;

  // hide  assignment operator
  G4VEnergyLossProcess(G4VEnergyLossProcess &) = delete;
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right) = delete;

  //------------------------------------------------------------------------
  // Get/set parameters to configure the process at initialisation time
  //------------------------------------------------------------------------

  // Add subcut processor for the region
  void ActivateSubCutoff(const G4Region* region);

  // Activate biasing
  void SetCrossSectionBiasingFactor(G4double f, G4bool flag = true);

  void ActivateForcedInteraction(G4double length, 
                                 const G4String& region,
                                 G4bool flag = true);

  void ActivateSecondaryBiasing(const G4String& region, G4double factor,
                                G4double energyLimit);

  inline void SetLossFluctuations(G4bool val);

  inline void SetSpline(G4bool val);
  inline void SetCrossSectionType(G4CrossSectionType val);
  inline G4CrossSectionType CrossSectionType() const;

  // Set/Get flag "isIonisation"
  void SetIonisation(G4bool val);
  inline G4bool IsIonisationProcess() const;

  // Redefine parameteters for stepping control
  void SetLinearLossLimit(G4double val);
  void SetStepFunction(G4double v1, G4double v2);
  void SetLowestEnergyLimit(G4double);

  inline G4int NumberOfSubCutoffRegions() const;

  //------------------------------------------------------------------------
  // Specific methods to path Physics Tables to the process
  //------------------------------------------------------------------------

  void SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType);
  void SetCSDARangeTable(G4PhysicsTable* pRange);
  void SetRangeTableForLoss(G4PhysicsTable* p);
  void SetInverseRangeTable(G4PhysicsTable* p);
  void SetLambdaTable(G4PhysicsTable* p);

  void SetTwoPeaksXS(std::vector<G4TwoPeaksXS*>*);
  void SetEnergyOfCrossSectionMax(std::vector<G4double>*);

  //------------------------------------------------------------------------
  // Specific methods to define custom Physics Tables to the process
  //------------------------------------------------------------------------

  // Binning for dEdx, range, inverse range and lambda tables
  void SetDEDXBinning(G4int nbins);

  // Min kinetic energy for tables
  void SetMinKinEnergy(G4double e);
  inline G4double MinKinEnergy() const;

  // Max kinetic energy for tables
  void SetMaxKinEnergy(G4double e);
  inline G4double MaxKinEnergy() const;

  // Biasing parameters
  inline G4double CrossSectionBiasingFactor() const;

  // Return values for given G4MaterialCutsCouple
  inline G4double GetDEDX(G4double kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetCSDADEDX(G4double kineticEnergy, 
                              const G4MaterialCutsCouple*);
  inline G4double GetDEDX(G4double kineticEnergy, const G4MaterialCutsCouple*,
                          G4double logKineticEnergy);
  inline G4double GetRange(G4double kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetRange(G4double kineticEnergy, const G4MaterialCutsCouple*,
                           G4double logKineticEnergy);
  inline G4double GetCSDARange(G4double kineticEnergy, 
                               const G4MaterialCutsCouple*);
  inline G4double GetKineticEnergy(G4double range, 
                                   const G4MaterialCutsCouple*);
  inline G4double GetLambda(G4double kineticEnergy,const G4MaterialCutsCouple*);
  inline G4double GetLambda(G4double kineticEnergy,const G4MaterialCutsCouple*,
                            G4double logKineticEnergy);

  inline G4bool TablesAreBuilt() const;

  // Access to specific tables
  inline G4PhysicsTable* DEDXTable() const;
  inline G4PhysicsTable* DEDXunRestrictedTable() const;
  inline G4PhysicsTable* IonisationTable() const;
  inline G4PhysicsTable* CSDARangeTable() const;
  inline G4PhysicsTable* RangeTableForLoss() const;
  inline G4PhysicsTable* InverseRangeTable() const;
  inline G4PhysicsTable* LambdaTable() const;
  inline std::vector<G4TwoPeaksXS*>* TwoPeaksXS() const;
  inline std::vector<G4double>* EnergyOfCrossSectionMax() const;

  inline G4bool UseBaseMaterial() const;

  //------------------------------------------------------------------------
  // Run time method for simulation of ionisation
  //------------------------------------------------------------------------

  // access atom on which interaction happens
  const G4Element* GetCurrentElement() const;

  // Set scaling parameters for ions is needed to G4EmCalculator
  void SetDynamicMassCharge(G4double massratio, G4double charge2ratio);

private:

  void FillSecondariesAlongStep(G4double weight);

  void PrintWarning(const G4String&, G4double val) const;

  // define material and indexes
  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Compute values using scaling relation, mass and charge of based particle
  //------------------------------------------------------------------------
  inline G4double GetDEDXForScaledEnergy(G4double scaledKinE);
  inline G4double GetDEDXForScaledEnergy(G4double scaledKinE,
                                         G4double logScaledKinE);
  inline G4double GetIonisationForScaledEnergy(G4double scaledKinE);
  inline G4double GetScaledRangeForScaledEnergy(G4double scaledKinE);
  inline G4double GetScaledRangeForScaledEnergy(G4double scaledKinE,
                                                G4double logScaledKinE);

  inline G4double GetLimitScaledRangeForScaledEnergy(G4double scaledKinE);
  inline G4double GetLimitScaledRangeForScaledEnergy(G4double scaledKinE,
                                                     G4double logScaledKinE);

  inline G4double ScaledKinEnergyForLoss(G4double range);
  inline G4double GetLambdaForScaledEnergy(G4double scaledKinE);
  inline G4double GetLambdaForScaledEnergy(G4double scaledKinE, 
                                           G4double logScaledKinE);

  void ComputeLambdaForScaledEnergy(G4double scaledKinE,
                                    G4double logScaledKinE);

  G4bool IsRegionForCubcutProcessor(const G4Track& aTrack);

protected:

  G4ParticleChangeForLoss     fParticleChange;
  const G4Material*           currentMaterial = nullptr;
  const G4MaterialCutsCouple* currentCouple = nullptr;

private:

  G4LossTableManager*         lManager;
  G4EmModelManager*           modelManager;
  G4VEmModel*                 currentModel = nullptr;
  G4EmBiasingManager*         biasManager = nullptr;
  G4SafetyHelper*             safetyHelper;
  G4EmParameters*             theParameters;  
  G4VEmFluctuationModel*      fluctModel = nullptr;
  G4VAtomDeexcitation*        atomDeexcitation = nullptr;
  G4VSubCutProducer*          subcutProducer = nullptr;

  const G4ParticleDefinition* particle = nullptr;
  const G4ParticleDefinition* baseParticle = nullptr;
  const G4ParticleDefinition* secondaryParticle = nullptr;
  G4EmDataHandler* theData = nullptr;

  G4PhysicsTable* theDEDXTable = nullptr;
  G4PhysicsTable* theDEDXunRestrictedTable = nullptr;
  G4PhysicsTable* theIonisationTable = nullptr;
  G4PhysicsTable* theRangeTableForLoss = nullptr;
  G4PhysicsTable* theCSDARangeTable = nullptr;
  G4PhysicsTable* theInverseRangeTable = nullptr;
  G4PhysicsTable* theLambdaTable = nullptr;

  std::vector<const G4Region*>* scoffRegions = nullptr;
  std::vector<G4VEmModel*>*     emModels = nullptr;
  const std::vector<G4int>*     theDensityIdx = nullptr;
  const std::vector<G4double>*  theDensityFactor = nullptr;
  const G4DataVector*           theCuts = nullptr;

  std::vector<G4double>* theEnergyOfCrossSectionMax = nullptr;
  std::vector<G4TwoPeaksXS*>* fXSpeaks = nullptr;

  G4double lowestKinEnergy;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;

  G4double linLossLimit = 0.01;
  G4double dRoverRange = 0.2;
  G4double finalRange;
  G4double lambdaFactor = 0.8;
  G4double invLambdaFactor;
  G4double biasFactor = 1.0;

  G4double massRatio = 1.0;
  G4double logMassRatio = 0.0;
  G4double fFactor = 1.0;
  G4double reduceFactor = 1.0;
  G4double chargeSqRatio = 1.0;
  G4double fRange = 0.0;
  G4double fRangeEnergy = 0.0;

protected:

  G4double preStepLambda = 0.0;
  G4double preStepKinEnergy = 0.0;
  G4double preStepLogKinEnergy = LOG_EKIN_MIN;
  G4double preStepScaledEnergy = 0.0;
  G4double preStepLogScaledEnergy = LOG_EKIN_MIN;
  G4double mfpKinEnergy = 0.0;

  std::size_t currentCoupleIndex = 0;

private:

  G4int nBins;
  G4int nBinsCSDA;
  G4int numberOfModels = 0;
  G4int nSCoffRegions = 0;
  G4int secID = _DeltaElectron;
  G4int tripletID = _TripletElectron;
  G4int biasID = _DeltaEBelowCut;
  G4int mainSecondaries = 1;

  std::size_t basedCoupleIndex = 0;
  std::size_t coupleIdxRange = 0;
  std::size_t idxDEDX = 0;
  std::size_t idxDEDXunRestricted = 0;
  std::size_t idxIonisation = 0;
  std::size_t idxRange = 0;
  std::size_t idxCSDA = 0;
  std::size_t idxSecRange = 0;
  std::size_t idxInverseRange = 0;
  std::size_t idxLambda = 0;

  G4GPILSelection aGPILSelection;
  G4CrossSectionType fXSType = fEmOnePeak;

  G4bool lossFluctuationFlag = true;
  G4bool rndmStepFlag = false;
  G4bool tablesAreBuilt = false;
  G4bool spline = true;
  G4bool isIon = false;
  G4bool isIonisation = true;
  G4bool useDeexcitation = false;
  G4bool biasFlag = false;
  G4bool weightFlag = false;
  G4bool isMaster = true;
  G4bool baseMat = false;
  G4bool actLinLossLimit = false;
  G4bool actLossFluc = false;
  G4bool actBinning = false;
  G4bool actMinKinEnergy = false;
  G4bool actMaxKinEnergy = false;

  std::vector<G4DynamicParticle*> secParticles;
  std::vector<G4Track*> scTracks;
};

// ======== Run time inline methods ================

inline std::size_t G4VEnergyLossProcess::CurrentMaterialCutsCoupleIndex() const 
{
  return currentCoupleIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SelectModel(G4double kinEnergy)
{
  currentModel = modelManager->SelectModel(kinEnergy, currentCoupleIndex);
  currentModel->SetCurrentCouple(currentCouple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::SelectModelForMaterial(
                   G4double kinEnergy, std::size_t& idx) const
{
  return modelManager->SelectModel(kinEnergy, idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEnergyLossProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple = couple;
    currentMaterial = couple->GetMaterial();
    basedCoupleIndex = currentCoupleIndex = couple->GetIndex();
    fFactor = chargeSqRatio*biasFactor;
    mfpKinEnergy = DBL_MAX;
    idxLambda = 0;
    if(baseMat) {
      basedCoupleIndex = (*theDensityIdx)[currentCoupleIndex];
      fFactor *= (*theDensityFactor)[currentCoupleIndex];
    }
    reduceFactor = 1.0/(fFactor*massRatio);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXForScaledEnergy(G4double e)
{
  /*
  G4cout << "G4VEnergyLossProcess::GetDEDX: Idx= " 
           << basedCoupleIndex << " E(MeV)= " << e 
         << " Emin= " << minKinEnergy << "  Factor= " << fFactor 
         << "  " << theDEDXTable << G4endl; */
  G4double x = fFactor*(*theDEDXTable)[basedCoupleIndex]->Value(e, idxDEDX);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

inline
G4double G4VEnergyLossProcess::GetDEDXForScaledEnergy(G4double e, G4double loge)
{
  /*
  G4cout << "G4VEnergyLossProcess::GetDEDX: Idx= " 
           << basedCoupleIndex << " E(MeV)= " << e 
         << " Emin= " << minKinEnergy << "  Factor= " << fFactor 
         << "  " << theDEDXTable << G4endl; */
  G4double x = fFactor*(*theDEDXTable)[basedCoupleIndex]->LogVectorValue(e,loge);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetIonisationForScaledEnergy(G4double e)
{
  G4double x = 
    fFactor*(*theIonisationTable)[basedCoupleIndex]->Value(e, idxIonisation);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetScaledRangeForScaledEnergy(G4double e)
{
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " lastIdx= " << lastIdx << "  " << theRangeTableForLoss << G4endl; 
  if(currentCoupleIndex != coupleIdxRange || fRangeEnergy != e) {
    coupleIdxRange = currentCoupleIndex;
    fRangeEnergy = e;
    fRange = reduceFactor*((*theRangeTableForLoss)[basedCoupleIndex])->Value(e, idxRange);
    if(e < minKinEnergy) { fRange *= std::sqrt(e/minKinEnergy); }
  }
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " R=  " << computedRange << "  " << theRangeTableForLoss << G4endl;
  return fRange;
}

inline G4double
G4VEnergyLossProcess::GetScaledRangeForScaledEnergy(G4double e, G4double loge)
{
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " lastIdx= " << lastIdx << "  " << theRangeTableForLoss << G4endl; 
  if(currentCoupleIndex != coupleIdxRange || fRangeEnergy != e) {
    coupleIdxRange = currentCoupleIndex;
    fRangeEnergy = e;
    fRange = reduceFactor*((*theRangeTableForLoss)[basedCoupleIndex])->LogVectorValue(e, loge);
    if(e < minKinEnergy) { fRange *= std::sqrt(e/minKinEnergy); }
  }
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " R=  " << fRange << "  " << theRangeTableForLoss << G4endl;
  return fRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLimitScaledRangeForScaledEnergy(G4double e)
{
  G4double x = ((*theCSDARangeTable)[basedCoupleIndex])->Value(e, idxCSDA);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLimitScaledRangeForScaledEnergy(G4double e,
                                                         G4double loge)
{
  G4double x = ((*theCSDARangeTable)[basedCoupleIndex])->LogVectorValue(e, loge);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::ScaledKinEnergyForLoss(G4double r)
{
  //G4cout << "G4VEnergyLossProcess::GetEnergy: Idx= " 
  //         << basedCoupleIndex << " R(mm)= " << r << "  " 
  //         << theInverseRangeTable << G4endl; 
  G4PhysicsVector* v = (*theInverseRangeTable)[basedCoupleIndex];
  G4double rmin = v->Energy(0);
  G4double e = 0.0; 
  if(r >= rmin) { e = v->Value(r, idxInverseRange); }
  else if(r > 0.0) {
    G4double x = r/rmin;
    e = minKinEnergy*x*x;
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambdaForScaledEnergy(G4double e)
{
  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e, idxLambda);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double
G4VEnergyLossProcess::GetLambdaForScaledEnergy(G4double e, G4double loge)
{
  return fFactor*((*theLambdaTable)[basedCoupleIndex])->LogVectorValue(e, loge);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetDEDX(G4double kinEnergy,
                              const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetDEDXForScaledEnergy(kinEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetDEDX(G4double kinEnergy,
                              const G4MaterialCutsCouple* couple,
                              G4double logKinEnergy)
{
  DefineMaterial(couple);
  return GetDEDXForScaledEnergy(kinEnergy*massRatio, logKinEnergy+logMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetRange(G4double kinEnergy,
                               const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetScaledRangeForScaledEnergy(kinEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetRange(G4double kinEnergy,
                               const G4MaterialCutsCouple* couple,
                               G4double logKinEnergy)
{
  DefineMaterial(couple);
  return GetScaledRangeForScaledEnergy(kinEnergy*massRatio, logKinEnergy+logMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetCSDARange(G4double kineticEnergy, 
                                   const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return (nullptr == theCSDARangeTable) ? DBL_MAX : 
    GetLimitScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetKineticEnergy(G4double range,
                                       const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return ScaledKinEnergyForLoss(range/reduceFactor)/massRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLambda(G4double kinEnergy,
                                const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return (nullptr != theLambdaTable) ? 
    GetLambdaForScaledEnergy(kinEnergy*massRatio) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLambda(G4double kinEnergy,
                                const G4MaterialCutsCouple* couple,
                                G4double logKinEnergy)
{
  DefineMaterial(couple);
  return (nullptr != theLambdaTable) ?
    GetLambdaForScaledEnergy(kinEnergy*massRatio, logKinEnergy+logMassRatio) 
    :  0.0;
}

// ======== Get/Set inline methods used at initialisation ================

inline void G4VEnergyLossProcess::SetFluctModel(G4VEmFluctuationModel* p)
{
  fluctModel = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmFluctuationModel* G4VEnergyLossProcess::FluctModel() const
{
  return fluctModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEnergyLossProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEnergyLossProcess::SetBaseParticle(const G4ParticleDefinition* p)
{
  baseParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEnergyLossProcess::Particle() const
{
  return particle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEnergyLossProcess::BaseParticle() const
{
  return baseParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* 
G4VEnergyLossProcess::SecondaryParticle() const
{
  return secondaryParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLossFluctuations(G4bool val)
{
  lossFluctuationFlag = val;
  actLossFluc = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetSpline(G4bool val)
{
  spline = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetCrossSectionType(G4CrossSectionType val)
{
  fXSType = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4CrossSectionType G4VEnergyLossProcess::CrossSectionType() const 
{
  return fXSType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossProcess::IsIonisationProcess() const
{
  return isIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossProcess::NumberOfSubCutoffRegions() const
{
  return nSCoffRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::CrossSectionBiasingFactor() const
{
  return biasFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossProcess::TablesAreBuilt() const
{
  return tablesAreBuilt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::DEDXTable() const
{
  return theDEDXTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::DEDXunRestrictedTable() const
{
  return theDEDXunRestrictedTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::IonisationTable() const
{
  return theIonisationTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::CSDARangeTable() const
{
  return theCSDARangeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::RangeTableForLoss() const
{
  return theRangeTableForLoss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::InverseRangeTable() const
{
  return theInverseRangeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossProcess::UseBaseMaterial() const
{
  return baseMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4double>* 
G4VEnergyLossProcess::EnergyOfCrossSectionMax() const
{
  return theEnergyOfCrossSectionMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4TwoPeaksXS*>* G4VEnergyLossProcess::TwoPeaksXS() const
{
  return fXSpeaks;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::size_t G4VEnergyLossProcess::NumberOfModels() const
{
  return numberOfModels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::EmModel(std::size_t index) const
{
  return (index < emModels->size()) ? (*emModels)[index] : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* 
G4VEnergyLossProcess::GetModelByIndex(std::size_t idx, G4bool ver) const
{
  return modelManager->GetModel((G4int)idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
