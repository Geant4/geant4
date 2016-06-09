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
// $Id: G4VEnergyLossProcess.hh,v 1.1 2003/11/12 16:18:09 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
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

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4VParticleChange;
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

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);
  // Build physics table during initialisation

  virtual void PrintInfoDefinition();

  // Print out of the class parameters

  G4PhysicsTable* BuildDEDXTable();

  G4PhysicsTable* BuildDEDXTableForPreciseRange();

  G4PhysicsTable* BuildLambdaTable();

  G4PhysicsTable* BuildLambdaSubTable();

  void SetParticles(const G4ParticleDefinition*,
                    const G4ParticleDefinition*);

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

  virtual void SetSubCutoff(G4bool) {};

  void SetDEDXTable(G4PhysicsTable* p);
  G4PhysicsTable* DEDXTable() const {return theDEDXTable;};

  void SetRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* RangeTable() const {return theRangeTable;};

  void SetInverseRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* InverseRangeTable() const {return theInverseRangeTable;};

  void SetSecondaryRangeTable(G4PhysicsTable* p);

  void SetLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* LambdaTable() {return theLambdaTable;};

  void SetSubLambdaTable(G4PhysicsTable* p);
  G4PhysicsTable* SubLambdaTable() {return theSubLambdaTable;};

  G4double GetDEDX(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetRange(G4double& kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetKineticEnergy(G4double& range, const G4MaterialCutsCouple* couple);

  G4double GetLambda(G4double kineticEnergy, const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process

  G4double GetDEDXDispersion(const G4MaterialCutsCouple *couple,
                             const G4DynamicParticle* dp,
                                   G4double& length);

  G4double MicroscopicCrossSection(G4double kineticEnergy,
                             const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process for a (energy, material)

  void SetLinearLossLimit(G4double val) {linLossLimit = val;};

  void SetLossFluctuations(G4bool val) {lossFluctuationFlag = val;};

  void SetIntegral(G4bool val);
  G4bool IsIntegral() const {return integral;}

  void SetRandomStep(G4bool val) {rndmStepFlag = val;};

  void SetMinSubRange(G4double val) {minSubRange = val;};

  void SetStepLimits(G4double v1, G4double v2);
  void SetStepFunction(G4double v1, G4double v2);

  G4bool TablesAreBuilt() const {return  tablesAreBuilt;};

  G4int NumberOfSubCutoffRegions() const {return nSCoffRegions;};

  G4double MeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimumStep,
                                     G4double& currentSafety);

  void      ResetNumberOfInteractionLengthLeft();
  // reset (determine the value of)NumberOfInteractionLengthLeft

protected:

  virtual
  G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual
  G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);

  virtual
  const G4ParticleDefinition* DefineBaseParticle(
          const G4ParticleDefinition*) {return 0;};

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

  G4VEmModel* SelectModel(G4double& kinEnergy);

  G4VSubCutoffProcessor* SubCutoffProcessor(size_t index);

  size_t CurrentMaterialCutsCoupleIndex() const {return currentMaterialIndex;};

  void SetMassRatio(G4double val) {massRatio = val;};

  void SetReduceFactor(G4double val) {reduceFactor = val;};

  void SetChargeSquare(G4double val) {chargeSquare = val;};

  void SetChargeSquareRatio(G4double val) {chargeSqRatio = val;};

private:

  void Clear();

  void DefineMaterial(const G4MaterialCutsCouple* couple);

  // hide  assignment operator

  G4VEnergyLossProcess(G4VEnergyLossProcess &);
  G4VEnergyLossProcess & operator=(const G4VEnergyLossProcess &right);

// =====================================================================

private:

  G4EmModelManager*                     modelManager;
  std::vector<G4VSubCutoffProcessor*>   scoffProcessors;
  std::vector<const G4Region*>          scoffRegions;
  G4int                                 nSCoffRegions;
  std::vector<G4int>                    idxSCoffRegions;

  // tables and vectors
  G4PhysicsTable*  theDEDXTable;
  G4PhysicsTable*  theRangeTable;
  G4PhysicsTable*  theSecondaryRangeTable;
  G4PhysicsTable*  theInverseRangeTable;
  G4PhysicsTable*  theLambdaTable;
  G4PhysicsTable*  theSubLambdaTable;
  G4double*        theDEDXAtMaxEnergy;
  G4double*        theRangeAtMaxEnergy;

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

  G4double faclow;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyForRange;
  G4double lowKinEnergy;
  G4double highKinEnergyForRange;

  G4double massRatio;
  G4double reduceFactor;
  G4double chargeSquare;
  G4double chargeSqRatio;

  G4double preStepLambda;
  G4double fRange;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double linLossLimit;
  G4double minSubRange;
  G4double dRoverRange;
  G4double finalRange;
  G4double defaultRoverRange;
  G4double defaultIntegralRange;

  G4bool lossFluctuationFlag;
  G4bool rndmStepFlag;
  G4bool hasRestProcess;
  G4bool tablesAreBuilt;
  G4bool integral;
  G4bool meanFreePath;
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
    if(integral && !meanFreePath) ResetNumberOfInteractionLengthLeft();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDX(G4double& kineticEnergy,
                                    const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  return ((*theDEDXTable)[currentMaterialIndex]->
          GetValue(kineticEnergy*massRatio, b))*chargeSqRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetRange(G4double& kineticEnergy,
                                            const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  G4double x;
  G4double e = kineticEnergy*massRatio;
  if (e < highKinEnergyForRange) {
    x = ((*theRangeTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = theRangeAtMaxEnergy[currentMaterialIndex] +
      (e - highKinEnergyForRange)/theDEDXAtMaxEnergy[currentMaterialIndex];
  }
  return x*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetKineticEnergy(G4double& range,
                                             const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  return ((*theInverseRangeTable)[currentMaterialIndex]->
         GetValue(range/reduceFactor, b))/massRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetDEDXDispersion(
                                  const G4MaterialCutsCouple *couple,
                                  const G4DynamicParticle* dp,
                                        G4double& length)
{
  DefineMaterial(couple);
  G4double tmax = MaxSecondaryEnergy(dp);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);
  return modelManager->GetDEDXDispersion(currentMaterial, dp, tmax, length,
                       currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetMeanFreePath(const G4Track& track,
                                                        G4double,
                                                        G4ForceCondition* cond)
{
  *cond = NotForced;

  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  if (meanFreePath) {
    G4bool b;
    preStepLambda = (((*theLambdaTable)[currentMaterialIndex])->
                      GetValue(preStepScaledEnergy, b)) * chargeSqRatio;
    if (integral) meanFreePath = false;
  }
  G4double x = DBL_MAX;
  if(0.0 < preStepLambda) x = 1.0/preStepLambda;
//  G4cout << GetProcessName() << ": e= " << preStepKinEnergy << " mfp= " << x << G4endl;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetContinuousStepLimit(const G4Track&,
                                  G4double, G4double currentMinStep, G4double&)
{
  G4double x = DBL_MAX;

  if (theRangeTable) {
    G4bool b;
    fRange = ((*theRangeTable)[currentMaterialIndex])->
            GetValue(preStepScaledEnergy, b)*reduceFactor;

    x = fRange;
    G4double y = x*dRoverRange;
    if(x > minStepLimit && y < currentMinStep ) {
      if (integral) x = std::max(y,minStepLimit);
      else {
        x = y + minStepLimit*(1.0 - dRoverRange)*(2.0 - minStepLimit/fRange);
        if(x > fRange) x = fRange;
        if(rndmStepFlag) x = minStepLimit + (x-minStepLimit)*G4UniformRand();
      }
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

inline G4VEmModel* G4VEnergyLossProcess::SelectModel(G4double& kinEnergy)
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

inline void G4VEnergyLossProcess::SetDEDXBinning(G4int nbins)
{
  nDEDXBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetDEDXBinningForPreciseRange(G4int nbins)
{
  nDEDXBinsForRange = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/*
inline G4int G4VEnergyLossProcess::DEDXBinning() const
{
  return nDEDXBins;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetLambdaBinning(G4int nbins)
{
  nLambdaBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/*
inline G4int G4VEnergyLossProcess::LambdaBinning() const
{
  return nLambdaBins;
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
  lowKinEnergy = minKinEnergy*faclow;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
  if(e < maxKinEnergyForRange) maxKinEnergyForRange = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossProcess::SetMaxKinEnergyForPreciseRange(G4double e)
{
  maxKinEnergyForRange = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossProcess::GetLambda(G4double kineticEnergy,
                                      const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  G4bool b;
  if(theLambdaTable) {
    G4double y = (((*theLambdaTable)[currentMaterialIndex])->
                    GetValue(kineticEnergy*massRatio, b));
    if(y > 0.0) x = 1.0/y;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VSubCutoffProcessor* G4VEnergyLossProcess::SubCutoffProcessor(size_t index)
{
  G4VSubCutoffProcessor* p = 0;
  if( nSCoffRegions ) p = scoffProcessors[idxSCoffRegions[index]];
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
