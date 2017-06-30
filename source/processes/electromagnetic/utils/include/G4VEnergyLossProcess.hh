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
// $Id: G4VEnergyLossProcess.hh 104349 2017-05-26 07:18:59Z gcosmo $
// GEANT4 tag $Name:
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
// Modifications:
//
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoffProcessors defined for regions (V.Ivanchenko)
// 17-02-03 Fix problem of store/restore tables (V.Ivanchenko)
// 26-02-03 Region dependent step limit (V.Ivanchenko)
// 26-03-03 Add GetDEDXDispersion (V.Ivanchenko)
// 09-04-03 Fix problem of negative range limit for non integral (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 21-07-03 Add UpdateEmModel method (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 14-01-04 Activate precise range calculation (V.Ivanchenko)
// 10-03-04 Fix problem of step limit calculation (V.Ivanchenko)
// 30-06-04 make destructor virtual (V.Ivanchenko)
// 05-07-04 fix problem of GenericIons seen at small cuts (V.Ivanchenko)
// 03-08-04 Add DEDX table to all processes for control on integral range(VI)
// 06-08-04 Clear up names of member functions (V.Ivanchenko)
// 27-08-04 Add NeedBuildTables method (V.Ivanchneko)
// 09-09-04 Bug fix for the integral mode with 2 peaks (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 11-04-05 Use MaxSecondaryEnergy from a model (V.Ivanchenko)
// 10-01-05 Remove SetStepLimits (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 13-01-06 Remove AddSubCutSecondaries and cleanup (V.Ivantchenko)
// 20-01-06 Introduce G4EmTableType and reducing number of methods (VI)
// 26-01-06 Add public method GetCSDARange (V.Ivanchenko)
// 22-03-06 Add SetDynamicMassCharge (V.Ivanchenko)
// 23-03-06 Use isIonisation flag (V.Ivanchenko)
// 13-05-06 Add method to access model by index (V.Ivanchenko)
// 14-01-07 add SetEmModel(index) and SetFluctModel() (mma)
// 15-01-07 Add separate ionisation tables and reorganise get/set methods for
//          dedx tables (V.Ivanchenko)
// 13-03-07 use SafetyHelper instead of navigator (V.Ivanchenko)
// 27-07-07 use stl vector for emModels instead of C-array (V.Ivanchenko)
// 25-09-07 More accurate handling zero xsect in 
//          PostStepGetPhysicalInteractionLength (V.Ivanchenko)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
// 15-07-08 Reorder class members for further multi-thread development (VI)
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
#include "G4UnitsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4EmTableType.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4EmParameters.hh"

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4Region;
class G4SafetyHelper;
class G4VAtomDeexcitation;
class G4VSubCutProducer;
class G4EmBiasingManager;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEnergyLossProcess : public G4VContinuousDiscreteProcess
{
public:

  G4VEnergyLossProcess(const G4String& name = "EnergyLoss",
                       G4ProcessType type = fElectromagnetic);

  virtual ~G4VEnergyLossProcess();

private:
  // clean vectors and arrays
  void Clean();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented in concrete processes
  //------------------------------------------------------------------------

public:
  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override = 0;
  
  virtual void PrintInfo() = 0;

  virtual void ProcessDescription(std::ostream& outFile) const; // = 0;

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut);

  //------------------------------------------------------------------------
  // Virtual methods implementation common to all EM ContinuousDiscrete 
  // processes. Further inheritance is not assumed 
  //------------------------------------------------------------------------

public:

  // prepare all tables
  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // build all tables
  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // build a table
  G4PhysicsTable* BuildDEDXTable(G4EmTableType tType = fRestricted);

  // build a table
  G4PhysicsTable* BuildLambdaTable(G4EmTableType tType = fRestricted);

  // summary printout after initialisation
  void PrintInfoDefinition(const G4ParticleDefinition& part);

  // Called before tracking of each new G4Track
  virtual void StartTracking(G4Track*) override;

  // Step limit from AlongStep 
  virtual G4double AlongStepGetPhysicalInteractionLength(
                                  const G4Track&,
                                  G4double  previousStepSize,
                                  G4double  currentMinimumStep,
                                  G4double& currentSafety,
                                  G4GPILSelection* selection) override;

  // Step limit from cross section
  virtual G4double PostStepGetPhysicalInteractionLength(
                                  const G4Track& track,
                                  G4double   previousStepSize,
                                  G4ForceCondition* condition) override;

  // AlongStep computations
  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, 
                                           const G4Step&) override;

  // Sampling of secondaries in vicinity of geometrical boundary
  // Return sum of secodaries energy 
  G4double SampleSubCutSecondaries(std::vector<G4Track*>&, const G4Step&, 
                                   G4VEmModel* model, G4int matIdx); 

  // PostStep sampling of secondaries
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, 
                                          const G4Step&) override;

  // Store all PhysicsTable in files.
  // Return false in case of any fatal failure at I/O  
  virtual G4bool StorePhysicsTable(const G4ParticleDefinition*,
                                   const G4String& directory,
                                   G4bool ascii = false) override;

  // Retrieve all Physics from a files.
  // Return true if all the Physics Table are built.
  // Return false if any fatal failure. 
  virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                                      const G4String& directory,
                                      G4bool ascii) override;

private:
  // store a table
  G4bool StoreTable(const G4ParticleDefinition* p, 
                    G4PhysicsTable*, G4bool ascii,
                    const G4String& directory, 
                    const G4String& tname);

  // retrieve a table
  G4bool RetrieveTable(const G4ParticleDefinition* p, 
                       G4PhysicsTable*, G4bool ascii,
                       const G4String& directory, 
                       const G4String& tname, 
                       G4bool mandatory);

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

  // access to cross section
  G4double MeanFreePath(const G4Track& track);

  // access to step limit
  G4double ContinuousStepLimit(const G4Track& track,
                               G4double previousStepSize,
                               G4double currentMinimumStep,
                               G4double& currentSafety);

protected:

  // implementation of the pure virtual method
  virtual G4double GetMeanFreePath(const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  // implementation of the pure virtual method
  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                          G4double previousStepSize,
                                          G4double currentMinimumStep,
                                          G4double& currentSafety) override;

  //------------------------------------------------------------------------
  // Run time method which may be also used by derived processes
  //------------------------------------------------------------------------

  // creeation of an empty vector for cross section
  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*, 
                                       G4double cut);

  inline size_t CurrentMaterialCutsCoupleIndex() const;

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  // Select model in run time
  inline void SelectModel(G4double kinEnergy);

public:
  // Select model by energy and region index
  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
                                            size_t& idx) const;

  // Add EM model coupled with fluctuation model for region, smaller value 
  // of order defines which pair of models will be selected for a given 
  // energy interval  
  void AddEmModel(G4int, G4VEmModel*, 
                  G4VEmFluctuationModel* fluc = 0,
                  const G4Region* region = nullptr);

  // Define new energy range for the model identified by the name
  void UpdateEmModel(const G4String&, G4double, G4double);

  // Assign a model to a process
  void SetEmModel(G4VEmModel*, G4int index=1);
  
  // return the assigned model
  G4VEmModel* EmModel(G4int index=1) const;
  
  // Access to models
  G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false) const;

  G4int NumberOfModels() const;

  // Assign a fluctuation model to a process
  void SetFluctModel(G4VEmFluctuationModel*);
  
  // return the assigned fluctuation model
  inline G4VEmFluctuationModel* FluctModel();
    
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

  //------------------------------------------------------------------------
  // Get/set parameters to configure the process at initialisation time
  //------------------------------------------------------------------------

  // Add subcutoff option for the region
  void ActivateSubCutoff(G4bool val, const G4Region* region = nullptr);

  // Activate biasing
  void SetCrossSectionBiasingFactor(G4double f, G4bool flag = true);

  void ActivateForcedInteraction(G4double length, 
                                 const G4String& region,
                                 G4bool flag = true);

  void ActivateSecondaryBiasing(const G4String& region, G4double factor,
                                G4double energyLimit);

  // Add subcutoff process (bremsstrahlung) to sample secondary 
  // particle production in vicinity of the geometry boundary
  void AddCollaborativeProcess(G4VEnergyLossProcess*);

  inline void SetLossFluctuations(G4bool val);

  inline void SetIntegral(G4bool val);
  inline G4bool IsIntegral() const;

  // Set/Get flag "isIonisation"
  void SetIonisation(G4bool val);
  inline G4bool IsIonisationProcess() const;

  // Redefine parameteters for stepping control
  void SetLinearLossLimit(G4double val);
  void SetStepFunction(G4double v1, G4double v2, G4bool lock=true);
  void SetLowestEnergyLimit(G4double);

  inline G4int NumberOfSubCutoffRegions() const;

  //------------------------------------------------------------------------
  // Specific methods to path Physics Tables to the process
  //------------------------------------------------------------------------

  void SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType);
  void SetCSDARangeTable(G4PhysicsTable* pRange);
  void SetRangeTableForLoss(G4PhysicsTable* p);
  void SetSecondaryRangeTable(G4PhysicsTable* p);
  void SetInverseRangeTable(G4PhysicsTable* p);

  void SetLambdaTable(G4PhysicsTable* p);
  void SetSubLambdaTable(G4PhysicsTable* p);

  // Binning for dEdx, range, inverse range and labda tables
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
  inline G4double GetDEDX(G4double& kineticEnergy, 
                          const G4MaterialCutsCouple*);
  inline G4double GetDEDXForSubsec(G4double& kineticEnergy, 
                                   const G4MaterialCutsCouple*);
  inline G4double GetRange(G4double& kineticEnergy, 
                           const G4MaterialCutsCouple*);
  inline G4double GetCSDARange(G4double& kineticEnergy, 
                               const G4MaterialCutsCouple*);
  inline G4double GetRangeForLoss(G4double& kineticEnergy, 
                                  const G4MaterialCutsCouple*);
  inline G4double GetKineticEnergy(G4double& range, 
                                   const G4MaterialCutsCouple*);
  inline G4double GetLambda(G4double& kineticEnergy, 
                            const G4MaterialCutsCouple*);

  inline G4bool TablesAreBuilt() const;

  // Access to specific tables
  inline G4PhysicsTable* DEDXTable() const;
  inline G4PhysicsTable* DEDXTableForSubsec() const;
  inline G4PhysicsTable* DEDXunRestrictedTable() const;
  inline G4PhysicsTable* IonisationTable() const;
  inline G4PhysicsTable* IonisationTableForSubsec() const;
  inline G4PhysicsTable* CSDARangeTable() const;
  inline G4PhysicsTable* SecondaryRangeTable() const;
  inline G4PhysicsTable* RangeTableForLoss() const;
  inline G4PhysicsTable* InverseRangeTable() const;
  inline G4PhysicsTable* LambdaTable() const;
  inline G4PhysicsTable* SubLambdaTable() const;

  //------------------------------------------------------------------------
  // Run time method for simulation of ionisation
  //------------------------------------------------------------------------

  // access atom on which interaction happens
  const G4Element* GetCurrentElement() const;

  // Set scaling parameters for ions is needed to G4EmCalculator
  inline void SetDynamicMassCharge(G4double massratio, G4double charge2ratio);

private:

  void FillSecondariesAlongStep(G4double& eloss, G4double& weight);

  void PrintWarning(G4String, G4double val);

  // define material and indexes
  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Compute values using scaling relation, mass and charge of based particle
  //------------------------------------------------------------------------

  inline G4double GetDEDXForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetSubDEDXForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetIonisationForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetSubIonisationForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetLimitScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  inline G4double ScaledKinEnergyForLoss(G4double range);
  inline G4double GetLambdaForScaledEnergy(G4double scaledKinEnergy);
  inline void ComputeLambdaForScaledEnergy(G4double scaledKinEnergy);

  // hide  assignment operator
  G4VEnergyLossProcess(G4VEnergyLossProcess &) = delete;
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right) = delete;

  // ======== Parameters of the class fixed at construction =========

  G4LossTableManager*         lManager;
  G4EmModelManager*           modelManager;
  G4EmBiasingManager*         biasManager;
  G4SafetyHelper*             safetyHelper;
  G4EmParameters*             theParameters;  

  const G4ParticleDefinition* secondaryParticle;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
  const G4ParticleDefinition* theGamma;
  const G4ParticleDefinition* theGenericIon;

  // ======== Parameters of the class fixed at initialisation =======

  std::vector<G4VEmModel*>              emModels;
  G4VEmFluctuationModel*                fluctModel;
  G4VAtomDeexcitation*                  atomDeexcitation;
  G4VSubCutProducer*                    subcutProducer;
  std::vector<const G4Region*>          scoffRegions;
  G4int                                 nSCoffRegions;
  G4bool*                               idxSCoffRegions;

  std::vector<G4VEnergyLossProcess*>    scProcesses;
  G4int                                 nProcesses;

  // tables and vectors
  G4PhysicsTable*             theDEDXTable;
  G4PhysicsTable*             theDEDXSubTable;
  G4PhysicsTable*             theDEDXunRestrictedTable;
  G4PhysicsTable*             theIonisationTable;
  G4PhysicsTable*             theIonisationSubTable;
  G4PhysicsTable*             theRangeTableForLoss;
  G4PhysicsTable*             theCSDARangeTable;
  G4PhysicsTable*             theSecondaryRangeTable;
  G4PhysicsTable*             theInverseRangeTable;
  G4PhysicsTable*             theLambdaTable;
  G4PhysicsTable*             theSubLambdaTable;

  size_t                      idxDEDX;
  size_t                      idxDEDXSub;
  size_t                      idxDEDXunRestricted;
  size_t                      idxIonisation;
  size_t                      idxIonisationSub;
  size_t                      idxRange;
  size_t                      idxCSDA;
  size_t                      idxSecRange;
  size_t                      idxInverseRange;
  size_t                      idxLambda;
  size_t                      idxSubLambda;

  std::vector<G4double>       theDEDXAtMaxEnergy;
  std::vector<G4double>       theRangeAtMaxEnergy;
  std::vector<G4double>       theEnergyOfCrossSectionMax;
  std::vector<G4double>       theCrossSectionMax;

  const std::vector<G4double>* theDensityFactor;
  const std::vector<G4int>*    theDensityIdx;

  const G4DataVector*         theCuts;
  const G4DataVector*         theSubCuts;

  const G4ParticleDefinition* baseParticle;

  G4int    nBins;
  G4int    nBinsCSDA;

  G4double lowestKinEnergy;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;

  G4double linLossLimit;
  G4double dRoverRange;
  G4double finalRange;
  G4double lambdaFactor;
  G4double biasFactor;

  G4bool   lossFluctuationFlag;
  G4bool   rndmStepFlag;
  G4bool   tablesAreBuilt;
  G4bool   integral;
  G4bool   isIon;
  G4bool   isIonisation;
  G4bool   useSubCutoff;
  G4bool   useDeexcitation;
  G4bool   biasFlag;
  G4bool   weightFlag;
  G4bool   isMaster;
  G4bool   actIntegral;
  G4bool   actStepFunc;
  G4bool   actLinLossLimit;
  G4bool   actLossFluc;
  G4bool   actBinning;
  G4bool   actMinKinEnergy;
  G4bool   actMaxKinEnergy;

protected:

  G4ParticleChangeForLoss          fParticleChange;

  // ======== Cached values - may be state dependent ================

private:

  std::vector<G4DynamicParticle*>  secParticles;
  std::vector<G4Track*>            scTracks;

  const G4ParticleDefinition* particle;

  G4VEmModel*                 currentModel;
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentCoupleIndex;
  size_t                      basedCoupleIndex;
  size_t                      lastIdx;

  G4double massRatio;
  G4double fFactor;
  G4double reduceFactor;
  G4double chargeSqRatio;

  G4double preStepLambda;
  G4double fRange;
  G4double computedRange;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double preStepRangeEnergy;
  G4double mfpKinEnergy;

  G4GPILSelection  aGPILSelection;

  G4int    secID;  
  G4int    subsecID;  
  G4int    biasID;  
};

// ======== Run time inline methods ================

inline size_t G4VEnergyLossProcess::CurrentMaterialCutsCoupleIndex() const 
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
                   G4double kinEnergy, size_t& idx) const
{
  return modelManager->SelectModel(kinEnergy, idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEnergyLossProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentCoupleIndex = couple->GetIndex();
    basedCoupleIndex   = (*theDensityIdx)[currentCoupleIndex];
    fFactor = chargeSqRatio*biasFactor*(*theDensityFactor)[currentCoupleIndex];
    reduceFactor = 1.0/(fFactor*massRatio);
    mfpKinEnergy = DBL_MAX;
    idxLambda = idxSubLambda = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetDynamicMassCharge(G4double massratio,
                                                       G4double charge2ratio)
{
  massRatio     = massratio;
  fFactor = charge2ratio*biasFactor*(*theDensityFactor)[currentCoupleIndex];
  chargeSqRatio = charge2ratio;
  reduceFactor  = 1.0/(fFactor*massRatio);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetSubDEDXForScaledEnergy(G4double e)
{
  G4double x = 
    fFactor*(*theDEDXSubTable)[basedCoupleIndex]->Value(e, idxDEDXSub);
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

inline 
G4double G4VEnergyLossProcess::GetSubIonisationForScaledEnergy(G4double e)
{
  G4double x = fFactor*
    (*theIonisationSubTable)[basedCoupleIndex]->Value(e, idxIonisationSub);
  if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetScaledRangeForScaledEnergy(G4double e)
{
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " lastIdx= " << lastIdx << "  " << theRangeTableForLoss << G4endl; 
  if(basedCoupleIndex != lastIdx || preStepRangeEnergy != e) {
    lastIdx = basedCoupleIndex;
    preStepRangeEnergy = e;
    computedRange = 
      ((*theRangeTableForLoss)[basedCoupleIndex])->Value(e, idxRange);
    if(e < minKinEnergy) { computedRange *= std::sqrt(e/minKinEnergy); }
  }
  //G4cout << "G4VEnergyLossProcess::GetScaledRange: Idx= " 
  //         << basedCoupleIndex << " E(MeV)= " << e 
  //         << " R=  " << fRange << "  " << theRangeTableForLoss << G4endl; 

  return computedRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLimitScaledRangeForScaledEnergy(G4double e)
{
  G4double x;
  if (e < maxKinEnergyCSDA) {
    x = ((*theCSDARangeTable)[basedCoupleIndex])->Value(e, idxCSDA);
    if(e < minKinEnergy) { x *= std::sqrt(e/minKinEnergy); }
  } else {
    x = theRangeAtMaxEnergy[basedCoupleIndex] +
      (e - maxKinEnergyCSDA)/theDEDXAtMaxEnergy[basedCoupleIndex];
  }
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
G4VEnergyLossProcess::GetDEDX(G4double& kineticEnergy,
                              const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetDEDXForScaledEnergy(kineticEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetDEDXForSubsec(G4double& kineticEnergy,
                                       const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetSubDEDXForScaledEnergy(kineticEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetRange(G4double& kineticEnergy,
                               const G4MaterialCutsCouple* couple)
{
  G4double x = fRange;
  DefineMaterial(couple);
  if(theCSDARangeTable) {
    x = GetLimitScaledRangeForScaledEnergy(kineticEnergy*massRatio)
      * reduceFactor;
  } else if(theRangeTableForLoss) {
    x = GetScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetCSDARange(G4double& kineticEnergy, 
                                   const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return (theCSDARangeTable) ? 
    GetLimitScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor
    : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetRangeForLoss(G4double& kineticEnergy,
                                      const G4MaterialCutsCouple* couple)
{
  //  G4cout << "GetRangeForLoss: Range from " << GetProcessName() << G4endl;
  DefineMaterial(couple);
  return GetScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetKineticEnergy(G4double& range,
                                       const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return ScaledKinEnergyForLoss(range/reduceFactor)/massRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4VEnergyLossProcess::GetLambda(G4double& kineticEnergy,
                                const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return theLambdaTable ? GetLambdaForScaledEnergy(kineticEnergy*massRatio) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::ComputeLambdaForScaledEnergy(G4double e)
{
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentCoupleIndex];
  if (e <= mfpKinEnergy) {
    preStepLambda = GetLambdaForScaledEnergy(e);

  } else {
    G4double e1 = e*lambdaFactor;
    if(e1 > mfpKinEnergy) {
      preStepLambda  = GetLambdaForScaledEnergy(e);
      G4double preStepLambda1 = GetLambdaForScaledEnergy(e1);
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

inline void G4VEnergyLossProcess::SetFluctModel(G4VEmFluctuationModel* p)
{
  fluctModel = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmFluctuationModel* G4VEnergyLossProcess::FluctModel()
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

inline void G4VEnergyLossProcess::SetIntegral(G4bool val)
{
  integral = val;
  actIntegral = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4bool G4VEnergyLossProcess::IsIntegral() const 
{
  return integral;
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
  return  tablesAreBuilt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::DEDXTable() const
{
  return theDEDXTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::DEDXTableForSubsec() const
{
  return theDEDXSubTable;
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

inline G4PhysicsTable* G4VEnergyLossProcess::IonisationTableForSubsec() const
{
  return theIonisationSubTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::CSDARangeTable() const
{
  return theCSDARangeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::SecondaryRangeTable() const
{
  return theSecondaryRangeTable;
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

inline G4PhysicsTable* G4VEnergyLossProcess::SubLambdaTable() const
{
  return theSubLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
