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
// $Id: G4VEnergyLossProcess.hh,v 1.57 2006/08/15 16:21:39 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
class G4Navigator;

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

  virtual std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double& tcut) = 0;

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------
protected:

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut);

  virtual void CorrectionsAlongStep(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
			           G4double& eloss,
                                   G4double& length);

  virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);

  //------------------------------------------------------------------------
  // Generic methods common to all processes 
  //------------------------------------------------------------------------
public:

  void PrintInfoDefinition();

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double SampleRange();

  G4PhysicsTable* BuildDEDXTable(G4EmTableType tType = fRestricted);

  G4PhysicsTable* BuildLambdaTable(G4EmTableType tType = fRestricted);

  void SetBaseParticle(const G4ParticleDefinition* p);

  const G4ParticleDefinition* Particle() const;
  const G4ParticleDefinition* BaseParticle() const;
  const G4ParticleDefinition* SecondaryParticle() const;

  // Binning for dEdx, range, and inverse range tables
  void SetDEDXBinning(G4int nbins);
  void SetLambdaBinning(G4int nbins);

  // Binning for dEdx, range, and inverse range tables
  void SetDEDXBinningForCSDARange(G4int nbins);

  // Min kinetic energy for tables
  void SetMinKinEnergy(G4double e);
  G4double MinKinEnergy() const;

  // Max kinetic energy for tables
  void SetMaxKinEnergy(G4double e);
  G4double MaxKinEnergy() const;

  // Max kinetic energy for tables
  void SetMaxKinEnergyForCSDARange(G4double e);

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

  // Add EM model coupled with fluctuation model for the region
  void AddEmModel(G4int, G4VEmModel*, G4VEmFluctuationModel* fluc = 0,
                                const G4Region* region = 0);

  // Define new energy range for the model identified by the name
  void UpdateEmModel(const G4String&, G4double, G4double);

  // Add subcutoff processor for the region
  void ActivateSubCutoff(G4bool val, const G4Region* region = 0);

  // Activate deexcitation code
  virtual void ActivateDeexcitation(G4bool, const G4Region* region = 0);

  void SetDEDXTable(G4PhysicsTable* p);
  G4PhysicsTable* DEDXTable() const;

  void SetDEDXTableForSubsec(G4PhysicsTable* p);
  G4PhysicsTable* DEDXTableForSubsec() const;

  void SetDEDXunRestrictedTable(G4PhysicsTable* p);
  G4PhysicsTable* DEDXunRestrictedTable() const;

  void SetCSDARangeTable(G4PhysicsTable* pRange);
  G4PhysicsTable* CSDARangeTable() const;

  void SetRangeTableForLoss(G4PhysicsTable* p);
  G4PhysicsTable* RangeTableForLoss() const;

  void SetInverseRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* InverseRangeTable() const;

  void SetSecondaryRangeTable(G4PhysicsTable* p);

  void SetLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* LambdaTable();

  void SetSubLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* SubLambdaTable();

  // Return values for given G4MaterialCutsCouple
  G4double GetDEDX(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  G4double GetDEDXForSubsec(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  G4double GetRange(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  G4double GetCSDARange(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  G4double GetRangeForLoss(G4double& kineticEnergy, const G4MaterialCutsCouple*);
  G4double GetKineticEnergy(G4double& range, const G4MaterialCutsCouple*);
  G4double GetLambda(G4double& kineticEnergy, const G4MaterialCutsCouple*);

  G4double GetDEDXDispersion(const G4MaterialCutsCouple *couple,
                             const G4DynamicParticle* dp,
                                   G4double length);

  G4double MicroscopicCrossSection(G4double kineticEnergy,
                             const G4MaterialCutsCouple* couple);

  void SetLossFluctuations(G4bool val);

  void SetRandomStep(G4bool val);

  void SetIntegral(G4bool val);
  G4bool IsIntegral() const;

  // Redefine parameteters for stepping control
  //
  void SetLinearLossLimit(G4double val);
  void SetMinSubRange(G4double val);
  void SetStepFunction(G4double v1, G4double v2);
  void SetLambdaFactor(G4double val);

  G4bool TablesAreBuilt() const;

  G4int NumberOfSubCutoffRegions() const;

  // Helper functions
  G4double MeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimumStep,
                                     G4double& currentSafety);

  // reset NumberOfInteractionLengthLeft
  void ResetNumberOfInteractionLengthLeft();

  G4VEmModel* SelectModelForMaterial(G4double kinEnergy, size_t& idx) const;

  // Set/Get flag "isIonisation"
  void SetIonisation(G4bool val);
  G4bool IsIonisationProcess() const;

  void AddCollaborativeProcess(G4VEnergyLossProcess*);

  void SampleSubCutSecondaries(std::vector<G4Track*>&, const G4Step&, 
                               G4double& cut, G4VEmModel* model); 

  // Set scaling parameters
  void SetDynamicMassCharge(G4double massratio, G4double charge2ratio);

  // Access to models
  G4VEmModel* GetModelByIndex(G4int idx = 0);

  G4int NumberOfModels();

protected:

  void SetParticle(const G4ParticleDefinition* p);

  void SetSecondaryParticle(const G4ParticleDefinition* p);

  G4VEmModel* SelectModel(G4double kinEnergy);

  size_t CurrentMaterialCutsCoupleIndex() const;

  G4double GetCurrentRange() const;

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*, G4double cut);

private:

  // Clear tables
  void Clear();

  void DefineMaterial(const G4MaterialCutsCouple* couple);

  // Returnd values for scaled energy and base particles mass
  //
  G4double GetDEDXForScaledEnergy(G4double scaledKinEnergy);
  G4double GetSubDEDXForScaledEnergy(G4double scaledKinEnergy);
  G4double GetScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  G4double GetLimitScaledRangeForScaledEnergy(G4double scaledKinEnergy);
  G4double GetLambdaForScaledEnergy(G4double scaledKinEnergy);
  G4double ScaledKinEnergyForLoss(G4double range);
  void ComputeLambdaForScaledEnergy(G4double scaledKinEnergy);

  // hide  assignment operator

  G4VEnergyLossProcess(G4VEnergyLossProcess &);
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right);

// =====================================================================

protected:

  G4ParticleChangeForLoss               fParticleChange;

private:

  G4EmModelManager*                     modelManager;
  std::vector<const G4Region*>          scoffRegions;
  G4int                                 nSCoffRegions;
  G4int*                                idxSCoffRegions;
  std::vector<G4Track*>                 scTracks;
  std::vector<G4VEnergyLossProcess*>    scProcesses;
  G4int                                 nProcesses;

  // tables and vectors
  G4PhysicsTable*             theDEDXTable;
  G4PhysicsTable*             theDEDXSubTable;
  G4PhysicsTable*             theDEDXunRestrictedTable;
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

  G4Navigator*                navigator;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* baseParticle;
  const G4ParticleDefinition* secondaryParticle;
  const G4ParticleDefinition* thePositron;

  G4PhysicsVector*            vstrag;

  // cash
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;

  G4int    nBins;
  G4int    nBinsCSDA;
  G4int    nWarnings;

  G4double lowestKinEnergy;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;

  G4double massRatio;
  G4double reduceFactor;
  G4double chargeSquare;
  G4double chargeSqRatio;

  G4double preStepLambda;
  G4double preStepMFP;
  G4double fRange;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double linLossLimit;
  G4double minSubRange;
  G4double dRoverRange;
  G4double finalRange;
  G4double lambdaFactor;
  G4double mfpKinEnergy;

  G4bool   lossFluctuationFlag;
  G4bool   lossFluctuationArePossible;
  G4bool   rndmStepFlag;
  G4bool   tablesAreBuilt;
  G4bool   integral;
  G4bool   meanFreePath;
  G4bool   aboveCSmax;
  G4bool   isIonisation;
  G4bool   useSubCutoff;
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
    if(!meanFreePath) ResetNumberOfInteractionLengthLeft();
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
  //  G4cout << "Range from " << GetProcessName() << "  e= " << kineticEnergy << " r= " << x << G4endl;
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
  G4double e = minKinEnergy;
  if(r <= rmin) {
    r /= rmin;
    e *= r*r;
  } else {
    G4bool b;
    e = v->GetValue(r, b);
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXDispersion(
                                  const G4MaterialCutsCouple *couple,
                                  const G4DynamicParticle* dp,
                                        G4double length)
{
  DefineMaterial(couple);
  G4double ekin = dp->GetKineticEnergy();
  G4VEmModel* currentModel = SelectModel(ekin*massRatio);
  G4double tmax = currentModel->MaxSecondaryKinEnergy(dp);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);
  G4double d = 0.0;
  G4VEmFluctuationModel* fm = currentModel->GetModelOfFluctuations();
  if(fm) d = fm->Dispersion(currentMaterial,dp,tmax,length);
  return d;
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
  meanFreePath  = false;
  aboveCSmax    = false;
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentMaterialIndex];
  if (e <= mfpKinEnergy) {
    preStepLambda = GetLambdaForScaledEnergy(e);
  } else {
    aboveCSmax  = true;
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

inline G4double G4VEnergyLossProcess::GetMeanFreePath(
                const G4Track& track, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  preStepKinEnergy = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  if(aboveCSmax && preStepScaledEnergy < mfpKinEnergy) 
    ResetNumberOfInteractionLengthLeft();
  DefineMaterial(track.GetMaterialCutsCouple());
  if (meanFreePath) {
    if (integral) ComputeLambdaForScaledEnergy(preStepScaledEnergy);
    else  preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy);
    if(0.0 < preStepLambda) preStepMFP = 1.0/preStepLambda;
    else                    preStepMFP = DBL_MAX;
  }
  //  G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy<< " eCSmax= " 
  //<<mfpKinEnergy<< " mfp= "<<preStepMFP<<G4endl;
  return preStepMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetContinuousStepLimit(const G4Track&,
                G4double, G4double currentMinStep, G4double&)
{
  G4double x = DBL_MAX;
  if(isIonisation) {
    fRange = GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;

    x = fRange;
    G4double y = x*dRoverRange;

    if(x > finalRange && y < currentMinStep ) {
      x = y + finalRange*(1.0 - dRoverRange)*(2.0 - finalRange/fRange);
      // G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
      // <<" range= "<<fRange <<" cMinSt="<<currentMinStep<< G4endl;
    } else if (rndmStepFlag) x = SampleRange();
  }
  //  G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
  //  <<" stepLimit= "<<x<<G4endl;
  return x;
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

inline void G4VEnergyLossProcess::ResetNumberOfInteractionLengthLeft()
{
  meanFreePath = true;
  aboveCSmax   = false;
  G4VProcess::ResetNumberOfInteractionLengthLeft();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MinPrimaryEnergy(
                const G4ParticleDefinition*, const G4Material*, G4double cut)
{
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::SelectModel(G4double kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::SelectModelForMaterial(
                   G4double kinEnergy, size_t& idx) const
{
  return modelManager->SelectModel(kinEnergy, idx);
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

inline void G4VEnergyLossProcess::CorrectionsAlongStep(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
			     G4double&,
			     G4double&)
{}

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
  chargeSquare  = charge2ratio*eplus*eplus;
  reduceFactor  = 1.0/(chargeSqRatio*massRatio);
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4double G4VEnergyLossProcess::GetCurrentRange() const
{
  return fRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4VEnergyLossProcess::AddCollaborativeProcess(
            G4VEnergyLossProcess* p)
{
  scProcesses.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::GetModelByIndex(G4int idx)
{
  return modelManager->GetModel(idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossProcess::NumberOfModels()
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
