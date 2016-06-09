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
// $Id: G4VEnergyLossProcess.hh,v 1.83 2008/09/12 16:19:01 vnivanch Exp $
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

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4Region;
class G4SafetyHelper;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEnergyLossProcess : public G4VContinuousDiscreteProcess
{
public:

  G4VEnergyLossProcess(const G4String& name = "EnergyLoss",
		       G4ProcessType type = fElectromagnetic);

  virtual ~G4VEnergyLossProcess();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented in concrete processes
  //------------------------------------------------------------------------

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;
  
  virtual void PrintInfo() = 0;

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

protected:

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
				    const G4Material*, G4double cut);

  //------------------------------------------------------------------------
  // Virtual methods common to all EM ContinuousDiscrete processes
  // Further inheritance is not assumed 
  //------------------------------------------------------------------------

public:

  void PrintInfoDefinition();

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
						 G4double  previousStepSize,
						 G4double  currentMinimumStep,
						 G4double& currentSafety,
						 G4GPILSelection* selection);

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
						G4double   previousStepSize,
						G4ForceCondition* condition);

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

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

protected:

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition);

  G4double GetContinuousStepLimit(const G4Track& track,
				  G4double previousStepSize,
				  G4double currentMinimumStep,
				  G4double& currentSafety);

  //------------------------------------------------------------------------
  // Specific methods for along/post step EM processes 
  //------------------------------------------------------------------------

public:

  void AddCollaborativeProcess(G4VEnergyLossProcess*);

  void SampleSubCutSecondaries(std::vector<G4Track*>&, const G4Step&, 
                               G4VEmModel* model, G4int matIdx,
			       G4double& extraEdep); 

  G4double GetDEDXDispersion(const G4MaterialCutsCouple *couple,
			     const G4DynamicParticle* dp,
			     G4double length);

  //------------------------------------------------------------------------
  // Specific methods to build and access Physics Tables
  //------------------------------------------------------------------------

  G4PhysicsTable* BuildDEDXTable(G4EmTableType tType = fRestricted);

  G4PhysicsTable* BuildLambdaTable(G4EmTableType tType = fRestricted);

  void SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType);
  void SetCSDARangeTable(G4PhysicsTable* pRange);
  void SetRangeTableForLoss(G4PhysicsTable* p);
  void SetInverseRangeTable(G4PhysicsTable* p);
  void SetSecondaryRangeTable(G4PhysicsTable* p);

  void SetLambdaTable(G4PhysicsTable* p);
  void SetSubLambdaTable(G4PhysicsTable* p);

  // Binning for dEdx, range, inverse range and labda tables
  inline void SetDEDXBinning(G4int nbins);
  inline void SetLambdaBinning(G4int nbins);

  // Binning for dEdx, range, and inverse range tables
  inline void SetDEDXBinningForCSDARange(G4int nbins);

  // Min kinetic energy for tables
  inline void SetMinKinEnergy(G4double e);
  inline G4double MinKinEnergy() const;

  // Max kinetic energy for tables
  inline void SetMaxKinEnergy(G4double e);
  inline G4double MaxKinEnergy() const;

  // Max kinetic energy for tables
  inline void SetMaxKinEnergyForCSDARange(G4double e);

  // Access to specific tables
  inline G4PhysicsTable* DEDXTable() const;
  inline G4PhysicsTable* DEDXTableForSubsec() const;
  inline G4PhysicsTable* DEDXunRestrictedTable() const;
  inline G4PhysicsTable* IonisationTable() const;
  inline G4PhysicsTable* IonisationTableForSubsec() const;
  inline G4PhysicsTable* CSDARangeTable() const;
  inline G4PhysicsTable* RangeTableForLoss() const;
  inline G4PhysicsTable* InverseRangeTable() const;
  inline G4PhysicsTable* LambdaTable();
  inline G4PhysicsTable* SubLambdaTable();

  // Return values for given G4MaterialCutsCouple
  inline G4double GetDEDX(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetDEDXForSubsec(G4double& kineticEnergy, 
				   const G4MaterialCutsCouple*);
  inline G4double GetRange(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetCSDARange(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetRangeForLoss(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  inline G4double GetKineticEnergy(G4double& range, const G4MaterialCutsCouple*);
  inline G4double GetLambda(G4double& kineticEnergy, const G4MaterialCutsCouple*);

  inline G4bool TablesAreBuilt() const;

  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

  inline void SetBaseParticle(const G4ParticleDefinition* p);
  inline const G4ParticleDefinition* Particle() const;
  inline const G4ParticleDefinition* BaseParticle() const;
  inline const G4ParticleDefinition* SecondaryParticle() const;

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  // Add EM model coupled with fluctuation model for the region
  inline void AddEmModel(G4int, G4VEmModel*, 
			 G4VEmFluctuationModel* fluc = 0,
			 const G4Region* region = 0);

  // Assign a model to a process
  inline void SetEmModel(G4VEmModel*, G4int index=1);
  
  // return the assigned model
  inline G4VEmModel* EmModel(G4int index=1);
  
  // Assign a fluctuation model to a process
  inline void SetFluctModel(G4VEmFluctuationModel*);
  
  // return the assigned fluctuation model
  inline G4VEmFluctuationModel* FluctModel();
    
  // Define new energy range for the model identified by the name
  inline void UpdateEmModel(const G4String&, G4double, G4double);

  // Access to models
  inline G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false);

  inline G4int NumberOfModels();

  //------------------------------------------------------------------------
  // Get/set parameters used for simulation of energy loss
  //------------------------------------------------------------------------

  inline void SetLossFluctuations(G4bool val);
  inline void SetRandomStep(G4bool val);
  inline void SetIntegral(G4bool val);
  inline G4bool IsIntegral() const;

  // Set/Get flag "isIonisation"
  inline void SetIonisation(G4bool val);
  inline G4bool IsIonisationProcess() const;

  // Redefine parameteters for stepping control
  //
  inline void SetLinearLossLimit(G4double val);
  inline void SetMinSubRange(G4double val);
  inline void SetStepFunction(G4double v1, G4double v2);
  inline void SetLambdaFactor(G4double val);


  // Add subcutoff option for the region
  void ActivateSubCutoff(G4bool val, const G4Region* region = 0);

  inline G4int NumberOfSubCutoffRegions() const;

  // Activate deexcitation code
  virtual void ActivateDeexcitation(G4bool, const G4Region* region = 0);

  //------------------------------------------------------------------------
  // Public interface to helper functions 
  //------------------------------------------------------------------------

  inline 
  G4VEmModel* SelectModelForMaterial(G4double kinEnergy, size_t& idx) const;

  inline G4double MeanFreePath(const G4Track& track);

  inline G4double ContinuousStepLimit(const G4Track& track,
				      G4double previousStepSize,
				      G4double currentMinimumStep,
				      G4double& currentSafety);

  //------------------------------------------------------------------------
  // Run time method for simulation of ionisation
  //------------------------------------------------------------------------

  // sample range at the end of a step
  inline G4double SampleRange();

  // Set scaling parameters for ions is needed to G4EmCalculator
  inline void SetDynamicMassCharge(G4double massratio, G4double charge2ratio);

  // Access to cross section table
  G4double CrossSectionPerVolume(G4double kineticEnergy,
				 const G4MaterialCutsCouple* couple);

protected:

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*, 
				       G4double cut);

  inline G4ParticleChangeForLoss* GetParticleChange();

  inline void SetParticle(const G4ParticleDefinition* p);

  inline void SetSecondaryParticle(const G4ParticleDefinition* p);

  inline void SelectModel(G4double kinEnergy);

  inline size_t CurrentMaterialCutsCoupleIndex() const;

  inline G4double GetCurrentRange() const;

private:

  //------------------------------------------------------------------------
  // Management of tables
  //------------------------------------------------------------------------

  void Clear();

  G4bool StoreTable(const G4ParticleDefinition* p, 
		    G4PhysicsTable*, G4bool ascii,
		    const G4String& directory, 
		    const G4String& tname);

  G4bool RetrieveTable(const G4ParticleDefinition* p, 
		       G4PhysicsTable*, G4bool ascii,
		       const G4String& directory, 
		       const G4String& tname, 
		       G4bool mandatory);

  // define material and indexes
  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  // Returnd values for scaled energy using mass of the base particle
  //
  inline G4double GetDEDXForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetSubDEDXForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetIonisationForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetSubIonisationForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetLimitScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  inline G4double GetLambdaForScaledEnergy(G4double scaledKinEnergy);
  inline G4double ScaledKinEnergyForLoss(G4double range);
  inline void ComputeLambdaForScaledEnergy(G4double scaledKinEnergy);

  // hide  assignment operator

  G4VEnergyLossProcess(G4VEnergyLossProcess &);
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right);

  // ======== Parameters of the class fixed at construction =========

  G4EmModelManager*           modelManager;
  G4SafetyHelper*             safetyHelper;

  const G4ParticleDefinition* secondaryParticle;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
  const G4ParticleDefinition* theGenericIon;

  G4PhysicsVector*            vstrag;

  // ======== Parameters of the class fixed at initialisation =======

  std::vector<G4VEmModel*>              emModels;
  G4VEmFluctuationModel*                fluctModel;
  std::vector<const G4Region*>          scoffRegions;
  G4int                                 nSCoffRegions;
  G4int*                                idxSCoffRegions;

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
  G4double*                   theDEDXAtMaxEnergy;
  G4double*                   theRangeAtMaxEnergy;
  G4double*                   theEnergyOfCrossSectionMax;
  G4double*                   theCrossSectionMax;

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
  G4double minSubRange;
  G4double dRoverRange;
  G4double finalRange;
  G4double lambdaFactor;

  G4bool   lossFluctuationFlag;
  G4bool   rndmStepFlag;
  G4bool   tablesAreBuilt;
  G4bool   integral;
  G4bool   isIon;
  G4bool   isIonisation;
  G4bool   useSubCutoff;

protected:

  G4ParticleChangeForLoss          fParticleChange;

  // ======== Cashed values - may be state dependent ================

private:

  std::vector<G4DynamicParticle*>  secParticles;
  std::vector<G4Track*>            scTracks;

  const G4ParticleDefinition* particle;

  G4VEmModel*                 currentModel;
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;

  G4int    nWarnings;

  G4double massRatio;
  G4double reduceFactor;
  G4double chargeSqRatio;

  G4double preStepLambda;
  G4double fRange;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double mfpKinEnergy;

  G4GPILSelection  aGPILSelection;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::DefineMaterial(
            const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    mfpKinEnergy = DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDX(G4double& kineticEnergy,
                                        const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetDEDXForScaledEnergy(kineticEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXForSubsec(G4double& kineticEnergy,
                                        const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetSubDEDXForScaledEnergy(kineticEnergy*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = 
    ((*theDEDXTable)[currentMaterialIndex]->GetValue(e, b))*chargeSqRatio;
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetSubDEDXForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = 
    ((*theDEDXSubTable)[currentMaterialIndex]->GetValue(e, b))*chargeSqRatio;
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetIonisationForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = 0.0;
  //  if(theIonisationTable) {
  x = ((*theIonisationTable)[currentMaterialIndex]->GetValue(e, b))
    *chargeSqRatio;
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  //}
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4VEnergyLossProcess::GetSubIonisationForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = 0.0;
  //if(theIonisationSubTable) {
  x = ((*theIonisationSubTable)[currentMaterialIndex]->GetValue(e, b))
    *chargeSqRatio;
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  //}
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRange(G4double& kineticEnergy,
                                         const G4MaterialCutsCouple* couple)
{
  G4double x = fRange;
  if(kineticEnergy != preStepKinEnergy || couple != currentCouple) { 
    DefineMaterial(couple);
    if(theCSDARangeTable)
      x = GetLimitScaledRangeForScaledEnergy(kineticEnergy*massRatio)
	* reduceFactor;
    else if(theRangeTableForLoss)
      x = GetScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetCSDARange(
       G4double& kineticEnergy, const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  if(theCSDARangeTable)
    x = GetLimitScaledRangeForScaledEnergy(kineticEnergy*massRatio)
      * reduceFactor;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLimitScaledRangeForScaledEnergy(
		G4double e)
{
  G4bool b;
  G4double x;

  if (e < maxKinEnergyCSDA) {
    x = ((*theCSDARangeTable)[currentMaterialIndex])->GetValue(e, b);
    if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  } else {
    x = theRangeAtMaxEnergy[currentMaterialIndex] +
         (e - maxKinEnergyCSDA)/theDEDXAtMaxEnergy[currentMaterialIndex];
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRangeForLoss(
                G4double& kineticEnergy,
		const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  if(theRangeTableForLoss) 
    x = GetScaledRangeForScaledEnergy(kineticEnergy*massRatio)*reduceFactor;
  //  G4cout << "Range from " << GetProcessName() 
  //         << "  e= " << kineticEnergy << " r= " << x << G4endl;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetScaledRangeForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = ((*theRangeTableForLoss)[currentMaterialIndex])->GetValue(e, b);
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetKineticEnergy(
                G4double& range,
		const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double r = range/reduceFactor;
  G4double e = ScaledKinEnergyForLoss(r)/massRatio;
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::ScaledKinEnergyForLoss(G4double r)
{
  G4PhysicsVector* v = (*theInverseRangeTable)[currentMaterialIndex];
  G4double rmin = v->GetLowEdgeEnergy(0);
  G4double e = 0.0; 
  if(r >= rmin) {
    G4bool b;
    e = v->GetValue(r, b);
  } else if(r > 0.0) {
    G4double x = r/rmin;
    e = minKinEnergy*x*x;
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambda(G4double& kineticEnergy,
					  const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = 0.0;
  if(theLambdaTable) x = GetLambdaForScaledEnergy(kineticEnergy*massRatio);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambdaForScaledEnergy(G4double e)
{
  G4bool b;
  return 
    chargeSqRatio*(((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::ComputeLambdaForScaledEnergy(G4double e)
{
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentMaterialIndex];
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
      preStepLambda = chargeSqRatio*theCrossSectionMax[currentMaterialIndex];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::ContinuousStepLimit(
         const G4Track& track, G4double x, G4double y, G4double& z)
{
  G4GPILSelection sel;
  return AlongStepGetPhysicalInteractionLength(track, x, y, z, &sel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::SampleRange()
{
  G4double e = amu_c2*preStepKinEnergy/particle->GetPDGMass();
  G4bool b;
  G4double s = fRange*std::pow(10.,vstrag->GetValue(e,b));
  G4double x = fRange + G4RandGauss::shoot(0.0,s);
  if(x > 0.0) fRange = x;
  return fRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MeanFreePath(const G4Track& track)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepLambda = GetLambdaForScaledEnergy(track.GetKineticEnergy()*massRatio);
  G4double x = DBL_MAX;
  if(DBL_MIN < preStepLambda) x = 1.0/preStepLambda;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MinPrimaryEnergy(
                const G4ParticleDefinition*, const G4Material*, G4double cut)
{
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SelectModel(G4double kinEnergy)
{
  currentModel = modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::SelectModelForMaterial(
                   G4double kinEnergy, size_t& idx) const
{
  return modelManager->SelectModel(kinEnergy, idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4ParticleChangeForLoss* G4VEnergyLossProcess::GetParticleChange()
{
  return &fParticleChange;
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

inline const G4ParticleDefinition* G4VEnergyLossProcess::SecondaryParticle() const
{
  return secondaryParticle;
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
  G4PhysicsTable* t = theDEDXTable;
  if(theIonisationTable) t = theIonisationTable; 
  return t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::IonisationTableForSubsec() const
{
  G4PhysicsTable* t = theDEDXSubTable;
  if(theIonisationSubTable) t = theIonisationSubTable; 
  return t;
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

inline G4PhysicsTable* G4VEnergyLossProcess::LambdaTable()
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEnergyLossProcess::SubLambdaTable()
{
  return theSubLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4bool G4VEnergyLossProcess::IsIntegral() const 
{
  return integral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline size_t G4VEnergyLossProcess::CurrentMaterialCutsCoupleIndex() const 
{
  return currentMaterialIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetDynamicMassCharge(G4double massratio,
                                                       G4double charge2ratio)
{
  massRatio     = massratio;
  chargeSqRatio = charge2ratio;
  reduceFactor  = 1.0/(chargeSqRatio*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4double G4VEnergyLossProcess::GetCurrentRange() const
{
  return fRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
void G4VEnergyLossProcess::AddEmModel(G4int order, G4VEmModel* p, 
				      G4VEmFluctuationModel* fluc,
				      const G4Region* region)
{
  modelManager->AddEmModel(order, p, fluc, region);
  if(p) p->SetParticleChange(pParticleChange, fluc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VEnergyLossProcess::GetModelByIndex(G4int idx, G4bool ver)
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossProcess::NumberOfModels()
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetEmModel(G4VEmModel* p, G4int index)
{
  G4int n = emModels.size();
  if(index >= n) for(G4int i=n; i<index+1; i++) {emModels.push_back(0);}
  emModels[index] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::EmModel(G4int index)
{
  G4VEmModel* p = 0;
  if(index >= 0 && index <  G4int(emModels.size())) p = emModels[index];
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

inline void G4VEnergyLossProcess::UpdateEmModel(const G4String& nam, 
						G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetIntegral(G4bool val)
{
  integral = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetBaseParticle(const G4ParticleDefinition* p)
{
  baseParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLinearLossLimit(G4double val)
{
  linLossLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLossFluctuations(G4bool val)
{
  lossFluctuationFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetRandomStep(G4bool val)
{
  rndmStepFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMinSubRange(G4double val)
{
  minSubRange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossProcess::TablesAreBuilt() const
{
  return  tablesAreBuilt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossProcess::NumberOfSubCutoffRegions() const
{
  return nSCoffRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetDEDXBinning(G4int nbins)
{
  nBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLambdaBinning(G4int nbins)
{
  nBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetDEDXBinningForCSDARange(G4int nbins)
{
  nBinsCSDA = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
  if(e < maxKinEnergyCSDA) maxKinEnergyCSDA = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMaxKinEnergyForCSDARange(G4double e)
{
  maxKinEnergyCSDA = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLambdaFactor(G4double val)
{
  if(val > 0.0 && val <= 1.0) lambdaFactor = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetIonisation(G4bool val)
{
  isIonisation = val;
  if(val) aGPILSelection = CandidateForSelection;
  else    aGPILSelection = NotCandidateForSelection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossProcess::IsIonisationProcess() const
{
  return isIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepFunction(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 0.999) dRoverRange = 1.0;
  currentCouple = 0;
  mfpKinEnergy  = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
