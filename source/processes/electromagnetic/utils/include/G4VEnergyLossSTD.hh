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
// File name:     G4BetheBlochModel
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

#ifndef G4VEnergyLossSTD_h
#define G4VEnergyLossSTD_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEffectiveChargeModel;
class G4VEmFluctuationModel;
class G4DataVector;
class G4VParticleChange;
class G4PhysicsTable;
class G4PhysicsVector;
class G4VSubCutoffProcessor;
class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEnergyLossSTD : public G4VContinuousDiscreteProcess
{
public:

  G4VEnergyLossSTD(const G4String& name = "EnergyLoss",
                         G4ProcessType type = fElectromagnetic);

 ~G4VEnergyLossSTD();

  void Initialise();

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
			           G4double& tmax,
			           G4double& eloss,
                                   G4double& kinEnergy) = 0;

  virtual void SecondariesPostStep(G4ParticleChange&,
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double& tcut,
                                   G4double& kinEnergy) = 0;

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);
    // True for all charged particles

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
  // Build physics table during initialisation

  virtual void PrintInfoDefinition() const;
  // Print out of the class parameters

  G4PhysicsTable* BuildDEDXTable();

  void SetParticles(const G4ParticleDefinition*,
                    const G4ParticleDefinition*,
                    const G4ParticleDefinition*);

  void SetParticle(const G4ParticleDefinition* p);
  void SetBaseParticle(const G4ParticleDefinition* p);
  void SetSecondaryParticle(const G4ParticleDefinition* p);

  const G4ParticleDefinition* Particle() const;
  const G4ParticleDefinition* BaseParticle() const;
  const G4ParticleDefinition* SecondaryParticle() const;
  // Print out of the class parameters

  void SetDEDXBinning(G4int nbins);
  G4int DEDXBinning() const;
    // Print out of the class parameters

  void SetLambdaBinning(G4int nbins);
  G4int LambdaBinning() const;
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

  void AddEmModel(G4VEmModel*, G4int, const G4Region* region = 0);

  void AddEmFluctuationModel(G4VEmFluctuationModel*, const G4Region* region = 0);

  virtual void SetSubCutoffProcessor(G4VSubCutoffProcessor*,
                               const G4Region* region = 0) {};

  virtual G4VSubCutoffProcessor* SubCutoffProcessor() {return 0;};

  virtual void SetSubCutoff(G4bool) {};

  void SetDEDXTable(G4PhysicsTable* p);
  G4PhysicsTable* DEDXTable() const {return theDEDXTable;};

  void SetRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* RangeTable() const {return theRangeTable;};

  void SetInverseRangeTable(G4PhysicsTable* p);
  G4PhysicsTable* InverseRangeTable() const {return theInverseRangeTable;};

  void SetSecondaryRangeTable(G4PhysicsTable* p);

  const G4PhysicsTable* LambdaTable() {return theLambdaTable;};

  G4double GetDEDX(G4double kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetRange(G4double kineticEnergy, const G4MaterialCutsCouple* couple);

  G4double GetKineticEnergy(G4double range, const G4MaterialCutsCouple* couple);

  G4double GetLambda(G4double kineticEnergy, const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process for a (energy, material)

  G4double MicroscopicCrossSection(G4double kineticEnergy,
                             const G4MaterialCutsCouple* couple);
  // It returns the MeanFreePath of the process for a (energy, material)

  void SetLinearLossLimit(G4double val) {linLossLimit = val;};

  void SetLossFluctuations(G4bool val) {lossFluctuationFlag = val;};

  void SetIntegral(G4bool val) {integral = val;};

  void SetRandomStep(G4bool val) {rndmStepFlag = val;};

  void SetMinSubRange(G4double val) {minSubRange = val;};

  void SetStepLimits(G4double v1, G4double v2);

  G4bool TablesAreBuilt() {return  tablesAreBuilt;};

  G4double MeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

  G4double ContinuousStepLimit(const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimumStep,
                                     G4double& currentSafety);

protected:

  virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);

  virtual const G4ParticleDefinition* DefineBaseParticle(
          const G4ParticleDefinition*) {return 0;};

  virtual G4PhysicsVector* DEDXPhysicsVector(const G4MaterialCutsCouple*);

  virtual G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  virtual G4PhysicsVector* SubLambdaPhysicsVector(const G4MaterialCutsCouple*);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut) = 0;

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp) = 0;

  G4VEmModel* SelectModel(G4double& kinEnergy);

  void SetMassRatio(G4double val) {massRatio = val;};

  void SetReduceFactor(G4double val) {reduceFactor = val;};

  void SetChargeSquare(G4double val) {chargeSquare = val;};

  void SetChargeSquareRatio(G4double val) {chargeSqRatio = val;};

private:

  void Clear();

  void DefineMaterial(const G4MaterialCutsCouple* couple);

  G4PhysicsTable* BuildLambdaTable();

  G4PhysicsTable* BuildLambdaSubTable();

  // hide  assignment operator

  G4VEnergyLossSTD(G4VEnergyLossSTD &);
  G4VEnergyLossSTD & operator=(const G4VEnergyLossSTD &right);

// =====================================================================

private:

  G4EmModelManager*                modelManager;
  G4VEmFluctuationModel*           emFluctModel;

  // tables and vectors
  G4PhysicsTable*  theDEDXTable;
  G4PhysicsTable*  theRangeTable;
  G4PhysicsTable*  theSecondaryRangeTable;
  G4PhysicsTable*  theInverseRangeTable;
  G4PhysicsTable*  theLambdaTable;

  const G4DataVector*    theCuts;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* baseParticle;
  const G4ParticleDefinition* secondaryParticle;
  const G4ParticleDefinition* theGamma;
  const G4ParticleDefinition* theElectron;

  // cash
  const G4Material*           currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;

  G4int    nDEDXBins;
  G4int    nLambdaBins;

  G4double minKinEnergy;
  G4double maxKinEnergy;

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
  //G4double maxFinalStep;
  G4double finalRange;
  G4double c1lim;
  G4double c2lim;
  G4double c3lim;

  G4bool lossFluctuationFlag;
  G4bool rndmStepFlag;
  G4bool hasRestProcess;
  G4bool tablesAreBuilt;
  G4bool integral;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEnergyLossSTD::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetDEDX(G4double kineticEnergy,
                                    const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  return ((*theDEDXTable)[currentMaterialIndex]->
          GetValue(kineticEnergy*massRatio, b))*chargeSqRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetRange(G4double kineticEnergy,
                                     const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  return ((*theRangeTable)[currentMaterialIndex]->
         GetValue(kineticEnergy*massRatio, b))*reduceFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetKineticEnergy(G4double range,
                                             const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4bool b;
  return ((*theInverseRangeTable)[currentMaterialIndex]->
         GetValue(range/reduceFactor, b))/massRatio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetMeanFreePath(const G4Track& track,
                                                        G4double,
                                                        G4ForceCondition* cond)
{
  G4bool b;
  *cond = NotForced;

  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy = track.GetKineticEnergy();
  preStepLambda = (((*theLambdaTable)[currentMaterialIndex])->
                      GetValue(preStepKinEnergy, b)) * chargeSqRatio;
  G4double x = DBL_MAX;
  if(0.0 < preStepLambda) x = 1.0/preStepLambda;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetContinuousStepLimit(const G4Track&,
                                               G4double, G4double, G4double&)
{
  G4double x = DBL_MAX;
  preStepScaledEnergy = preStepKinEnergy*massRatio;

  if(theRangeTable) {
    G4bool b;
    fRange = ((*theRangeTable)[currentMaterialIndex])->
                GetValue(preStepScaledEnergy, b)*reduceFactor;

    if(integral || fRange <= finalRange) {
      x = fRange;

    } else {
      x = c1lim*fRange+c2lim+c3lim/fRange;
      if(rndmStepFlag) x = finalRange + (x-finalRange)*G4UniformRand();
      if(x > fRange) x = fRange;
    }

  }

  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEnergyLossSTD::SelectModel(G4double& kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  if(!baseParticle) baseParticle = DefineBaseParticle(particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetBaseParticle(const G4ParticleDefinition* p)
{
  baseParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEnergyLossSTD::Particle() const
{
  return particle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEnergyLossSTD::BaseParticle() const
{
  return baseParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VEnergyLossSTD::SecondaryParticle() const
{
  return secondaryParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetDEDXBinning(G4int nbins)
{
  nDEDXBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossSTD::DEDXBinning() const
{
  return nDEDXBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetLambdaBinning(G4int nbins)
{
  nLambdaBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEnergyLossSTD::LambdaBinning() const
{
  return nLambdaBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEnergyLossSTD::GetLambda(G4double kineticEnergy,
                                      const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4double x = DBL_MAX;
  G4bool b;
  if(theLambdaTable) {
    G4double y = (((*theLambdaTable)[currentMaterialIndex])->
                    GetValue(kineticEnergy, b));
    if(y > 0.0) x = 1.0/y;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEnergyLossSTD::SetStepLimits(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  c1lim=dRoverRange;
  c2lim=2.*(1-dRoverRange)*finalRange;
  c3lim=-(1.-dRoverRange)*finalRange*finalRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
