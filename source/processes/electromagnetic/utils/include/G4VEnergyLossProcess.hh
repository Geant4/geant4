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
// $Id: G4VEnergyLossProcess.hh,v 1.21 2004/05/17 09:46:56 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02 $
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

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4PhysicsTable;
class G4PhysicsVector;
class G4VSubCutoffProcessor;
class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEnergyLossProcess : public G4VContinuousDiscreteProcess
{
public:

  G4VEnergyLossProcess(const G4String& name = "EnergyLoss",
                         G4ProcessType type = fElectromagnetic);

 ~G4VEnergyLossProcess();

  void Initialise();

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
			           G4double& tmax,
			           G4double& eloss,
                                   G4double& kinEnergy) = 0;

  virtual void SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double& tcut,
                                   G4double& kinEnergy) = 0;

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;
    // True for all charged particles

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
  // Build physics table during initialisation

  virtual void PrintInfoDefinition();

  // Print out of the class parameters

  G4PhysicsTable* BuildDEDXTable();

  G4PhysicsTable* BuildDEDXTableForPreciseRange();

  G4PhysicsTable* BuildLambdaTable();

  G4PhysicsTable* BuildLambdaSubTable();

  void SetParticle(const G4ParticleDefinition* p);
  void SetBaseParticle(const G4ParticleDefinition* p);
  void SetSecondaryParticle(const G4ParticleDefinition* p);

  const G4ParticleDefinition* Particle() const;
  const G4ParticleDefinition* BaseParticle() const;
  const G4ParticleDefinition* SecondaryParticle() const;
  // Particle definition

  void SetDEDXBinning(G4int nbins);
  // Binning for dEdx, range, and inverse range tables

  void SetDEDXBinningForPreciseRange(G4int nbins);
  // Binning for dEdx, range, and inverse range tables

  void SetLambdaBinning(G4int nbins);
  // Binning for lambda table

  void SetMinKinEnergy(G4double e);
  G4double MinKinEnergy() const;
  // Min kinetic energy for tables

  void SetMaxKinEnergy(G4double e);
  G4double MaxKinEnergy() const;
  // Max kinetic energy for tables

  void SetMaxKinEnergyForPreciseRange(G4double e);
  // Max kinetic energy for tables

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

  void AddEmModel(G4int, G4VEmModel*, G4VEmFluctuationModel* fluc = 0,
                                const G4Region* region = 0);
  // Add EM model coupled with fluctuation model for the region

  void UpdateEmModel(const G4String&, G4double, G4double);
  // Define new energy range for thhe model identified by the name

  void AddSubCutoffProcessor(G4VSubCutoffProcessor*, const G4Region* region = 0);
  // Add subcutoff processor for the region

  virtual void ActivateFluorescence(G4bool, const G4Region* region = 0);
  virtual void ActivateAugerElectronProduction(G4bool, const G4Region* region = 0);
  // Activate deexcitation code

  virtual void SetSubCutoff(G4bool);

  void SetDEDXTable(G4PhysicsTable* p);
  G4PhysicsTable* DEDXTable() const;

  void SetPreciseRangeTable(G4PhysicsTable* pRange);
  G4PhysicsTable* PreciseRangeTable() const;

  void SetRangeTableForLoss(G4PhysicsTable* p);
  G4PhysicsTable* RangeTableForLoss() const;

  void SetInverseRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* InverseRangeTable() const;

  void SetSecondaryRangeTable(G4PhysicsTable* p);

  void SetLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* LambdaTable();

  void SetSubLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* SubLambdaTable();

  G4double GetDEDX(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetRange(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetRangeForLoss(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetKineticEnergy(G4double& range, const G4MaterialCutsCouple* couple);

  G4double GetLambda(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process

  G4double GetDEDXDispersion(const G4MaterialCutsCouple *couple,
                             const G4DynamicParticle* dp,
                                   G4double length);

  G4double MicroscopicCrossSection(G4double kineticEnergy,
                             const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process for a (energy, material)

  void SetLinearLossLimit(G4double val);

  void SetLossFluctuations(G4bool val);

  void SetIntegral(G4bool val);
  G4bool IsIntegral() const;

  void SetRandomStep(G4bool val);

  void SetMinSubRange(G4double val);

  void SetStepLimits(G4double v1, G4double v2);

  void SetStepFunction(G4double v1, G4double v2);

  void SetLambdaFactor(G4double val);

  G4bool TablesAreBuilt() const;

  G4int NumberOfSubCutoffRegions() const;

  G4double MeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimumStep,
                                     G4double& currentSafety);

  void ResetNumberOfInteractionLengthLeft();
  // reset (determine the value of)NumberOfInteractionLengthLeft

protected:

  virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);

  virtual
  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition*);

  virtual
  G4PhysicsVector* DEDXPhysicsVector(const G4MaterialCutsCouple*);

  virtual
  G4PhysicsVector* DEDXPhysicsVectorForPreciseRange(const G4MaterialCutsCouple*);

  virtual
  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  virtual
  G4PhysicsVector* SubLambdaPhysicsVector(const G4MaterialCutsCouple*);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut) = 0;

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp) = 0;

  G4VEmModel* SelectModel(G4double kinEnergy);

  G4VSubCutoffProcessor* SubCutoffProcessor(size_t index);

  size_t CurrentMaterialCutsCoupleIndex() const;

  void SetMassRatio(G4double val);

  void SetReduceFactor(G4double val);

  void SetChargeSquare(G4double val);

  void SetChargeSquareRatio(G4double val);
  
  G4double GetCurrentRange() const;

private:

  void Clear();

  void DefineMaterial(const G4MaterialCutsCouple* couple);

  G4double GetDEDXForLoss(G4double kineticEnergy);

  G4double GetRangeForLoss(G4double kineticEnergy);

  G4double GetPreciseRange(G4double kineticEnergy);

  G4double GetLambda(G4double scaledKinEnergy);

  void ComputeLambda(G4double scaledKinEnergy);

  G4double ScaledKinEnergyForLoss(G4double range);

  // hide  assignment operator

  G4VEnergyLossProcess(G4VEnergyLossProcess &);
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right);

// =====================================================================

protected:

  G4ParticleChangeForLoss               fParticleChange;

private:

  G4EmModelManager*                     modelManager;
  std::vector<G4VSubCutoffProcessor*>   scoffProcessors;
  std::vector<const G4Region*>          scoffRegions;
  G4int                                 nSCoffRegions;
  std::vector<G4int>                    idxSCoffRegions;

  // tables and vectors
  G4PhysicsTable*             theDEDXTable;
  G4PhysicsTable*             theRangeTableForLoss;
  G4PhysicsTable*             thePreciseRangeTable;
  G4PhysicsTable*             theSecondaryRangeTable;
  G4PhysicsTable*             theInverseRangeTable;
  G4PhysicsTable*             theLambdaTable;
  G4PhysicsTable*             theSubLambdaTable;
  G4double*                   theDEDXAtMaxEnergy;
  G4double*                   theRangeAtMaxEnergy;
  G4double*                   theEnergyOfCrossSectionMax;
  G4double*                   theCrossSectionMax;

  const G4DataVector*         theCuts;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* baseParticle;
  const G4ParticleDefinition* secondaryParticle;

  // cash
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;
  G4double                    minStepLimit;

  G4int    nDEDXBins;
  G4int    nDEDXBinsForRange;
  G4int    nLambdaBins;

  G4double lowestKinEnergy;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyForRange;

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
  G4double defaultRoverRange;
  G4double defaultIntegralRange;
  G4double lambdaFactor;
  G4double mfpKinEnergy;

  G4bool   lossFluctuationFlag;
  G4bool   rndmStepFlag;
  G4bool   hasRestProcess;
  G4bool   tablesAreBuilt;
  G4bool   integral;
  G4bool   meanFreePath;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    minStepLimit = std::min(finalRange,
           currentCouple->GetProductionCuts()->GetProductionCut(idxG4ElectronCut));
    if(integral && (!meanFreePath || preStepScaledEnergy < mfpKinEnergy))
      ResetNumberOfInteractionLengthLeft();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDX(G4double& kineticEnergy,
                                        const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetDEDXForLoss(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXForLoss(G4double e)
{
  G4bool b;
  e *= massRatio;
  G4double x = ((*theDEDXTable)[currentMaterialIndex]->GetValue(e, b))*chargeSqRatio;
  if(e < minKinEnergy) x *= sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRange(G4double& kineticEnergy,
                                         const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  if(thePreciseRangeTable)       x = GetPreciseRange(kineticEnergy);
  else if(theRangeTableForLoss)  x = GetRangeForLoss(kineticEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetPreciseRange(G4double e)
{
  G4bool b;
  G4double x;
  e *= massRatio;

  if (e < maxKinEnergyForRange) {
    x = ((*thePreciseRangeTable)[currentMaterialIndex])->GetValue(e, b);
    if(e < minKinEnergy) x *= sqrt(e/minKinEnergy);

  } else {
    x = theRangeAtMaxEnergy[currentMaterialIndex] +
         (e - maxKinEnergyForRange)/theDEDXAtMaxEnergy[currentMaterialIndex];
  }
  return x*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRangeForLoss(G4double& kineticEnergy,
                                                const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  if(theRangeTableForLoss) x = GetRangeForLoss(kineticEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRangeForLoss(G4double e)
{
  G4bool b;
  e *= massRatio;
  G4double x = ((*theRangeTableForLoss)[currentMaterialIndex])->GetValue(e, b);
  if(e < minKinEnergy) x *= sqrt(e/minKinEnergy);
  return x*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetKineticEnergy(G4double& range,
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
  G4double tmax = MaxSecondaryEnergy(dp);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);
  return modelManager->GetDEDXDispersion(currentMaterial, dp, tmax, length,
                       currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambda(G4double& kineticEnergy,
                                          const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = 0.0;
  if(theLambdaTable) x = GetLambda(kineticEnergy*massRatio);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambda(G4double e)
{
  G4bool b;
  return chargeSqRatio*(((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::ComputeLambda(G4double e)
{
  meanFreePath  = false;
  mfpKinEnergy  = 0.0;
  G4double emax = theEnergyOfCrossSectionMax[currentMaterialIndex];
  if (e <= emax) preStepLambda = GetLambda(e);
  else {
    e *= lambdaFactor;
    if(e > emax) {
      mfpKinEnergy = e;
      preStepLambda = GetLambda(e);
    } else preStepLambda = theCrossSectionMax[currentMaterialIndex];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetMeanFreePath(const G4Track& track,
                                                            G4double, G4ForceCondition*)
{
  preStepKinEnergy = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  DefineMaterial(track.GetMaterialCutsCouple());
  if (meanFreePath) {
    if (integral) ComputeLambda(preStepScaledEnergy);
    else          preStepLambda = GetLambda(preStepScaledEnergy);
    if(0.0 < preStepLambda) preStepMFP = 1.0/preStepLambda;
    else                    preStepMFP = DBL_MAX;
  }
  // G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy<<" mfp= "<<preStepMFP<<G4endl;
  return preStepMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetContinuousStepLimit(const G4Track&,
                                      G4double, G4double currentMinStep, G4double&)
{
  G4double x = DBL_MAX;
  if(theRangeTableForLoss) {
    fRange = GetRange(preStepKinEnergy, currentCouple);

    x = fRange;
    G4double y = x*dRoverRange;

    if(x > minStepLimit && y < currentMinStep ) {

      x = y + minStepLimit*(1.0 - dRoverRange)*(2.0 - minStepLimit/fRange);
      //if(x >fRange || x<minStepLimit) G4cout << "!!! StepLimit problem!!!" << G4endl;
      //if(rndmStepFlag) x = minStepLimit + G4UniformRand()*(x-minStepLimit);
    }
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::ResetNumberOfInteractionLengthLeft()
{
  meanFreePath = true;
  G4VProcess::ResetNumberOfInteractionLengthLeft();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossProcess::SelectModel(G4double kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
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

inline G4VSubCutoffProcessor* G4VEnergyLossProcess::SubCutoffProcessor(size_t index)
{
  G4VSubCutoffProcessor* p = 0;
  if( nSCoffRegions ) p = scoffProcessors[idxSCoffRegions[index]];
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4PhysicsTable* G4VEnergyLossProcess::DEDXTable() const 
{
  return theDEDXTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4PhysicsTable* G4VEnergyLossProcess::PreciseRangeTable() const
{
  return thePreciseRangeTable;
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
  
inline void G4VEnergyLossProcess::SetMassRatio(G4double val) 
{
  massRatio = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4VEnergyLossProcess::SetReduceFactor(G4double val) 
{
  reduceFactor = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4VEnergyLossProcess::SetChargeSquare(G4double val) 
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4VEnergyLossProcess::SetChargeSquareRatio(G4double val) 
{
  chargeSqRatio = val;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline G4double G4VEnergyLossProcess::GetCurrentRange() const 
{
  return fRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
