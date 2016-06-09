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
// $Id: G4VEmProcess.hh,v 1.28 2005/05/12 11:06:52 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
// 16-09-04 Add flag for LambdaTable and method RecalculateLambda (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
// 09-05-05 Fix problem in logic when path boundary between materials (V.Ivantchenko)
//
// Class Description:
//
// It is the unified Discrete process

// -------------------------------------------------------------------
//

#ifndef G4VEmProcess_h
#define G4VEmProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForGamma.hh"

class G4Step;
class G4VEmModel;
class G4DataVector;
class G4VParticleChange;
class G4PhysicsTable;
class G4PhysicsVector;

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

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;

  virtual void PrintInfo() = 0;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

  virtual std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                 G4VEmModel*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  virtual G4double RecalculateLambda(G4double kinEnergy,
                               const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Generic methods common to all processes 
  //------------------------------------------------------------------------

public:

  void PrintInfoDefinition();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  void PreparePhysicsTable(const G4ParticleDefinition&);
  // Initialise for build of tables

  void BuildPhysicsTable(const G4ParticleDefinition&);
  // Build physics table during initialisation

  void SetLambdaBinning(G4int nbins);
  G4int LambdaBinning() const;
  // Binning for lambda table

  void SetMinKinEnergy(G4double e);
  G4double MinKinEnergy() const;
  // Min kinetic energy for tables

  void SetMaxKinEnergy(G4double e);
  G4double MaxKinEnergy() const;
  // Max kinetic energy for tables

  void SetLambdaFactor(G4double val);

  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                                 G4bool ascii = false);
    // Store PhysicsTable in a file.
    // Return false in case of failure at I/O

  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                                    G4bool ascii);
    // Retrieve Physics from a file.
    // (return true if the Physics Table can be build by using file)
    // (return false if the process has no functionality or in case of failure)
    // File name should is constructed as processName+particleName and the
    // should be placed under the directory specifed by the argument.

  void AddEmModel(G4int, G4VEmModel*, const G4Region* region = 0);
  // Add EM model coupled for the region

  void UpdateEmModel(const G4String&, G4double, G4double);
  // Define new energy range for the model identified by the name

  G4double GetLambda(G4double& kinEnergy, const G4MaterialCutsCouple* couple);
  // It returns the Lambda of the process

  const G4PhysicsTable* LambdaTable() const;

  G4double MicroscopicCrossSection(G4double kineticEnergy,
                             const G4MaterialCutsCouple* couple);
  // It returns the cross section of the process for energy/ material

  G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, G4double Z);
  // It returns the cross section of the process per atom

  G4double MeanFreePath(     const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition);

  const G4ParticleDefinition* Particle() const;
  const G4ParticleDefinition* SecondaryParticle() const;

  void ActivateDeexcitation(G4bool, const G4Region* r = 0);

  G4VEmModel* SelectModelForMaterial(G4double kinEnergy, size_t& idxRegion) const;

  void SetIntegral(G4bool val);
  G4bool IsIntegral() const;

  void SetApplyCuts(G4bool val);
  
protected:

  void SetParticle(const G4ParticleDefinition* p);
  
  void SetSecondaryParticle(const G4ParticleDefinition* p);

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition);

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  G4VEmModel* SelectModel(G4double& kinEnergy);

  size_t CurrentMaterialCutsCoupleIndex() const {return currentMaterialIndex;};

  void ResetNumberOfInteractionLengthLeft();

  G4double GetGammaEnergyCut();
  G4double GetElectronEnergyCut();

  void SetBuildTableFlag(G4bool val);

  void SetStartFromNullFlag(G4bool val);

private:

  void Clear();

  void DefineMaterial(const G4MaterialCutsCouple* couple);

  void ComputeIntegralLambda(G4double kinEnergy);

  G4double GetLambdaFromTable(G4double kinEnergy);

  G4double GetCurrentLambda(G4double kinEnergy);

  G4double ComputeCurrentLambda(G4double kinEnergy);

  void BuildLambdaTable();

  void FindLambdaMax();

  // hide  assignment operator

  G4VEmProcess(G4VEmProcess &);
  G4VEmProcess & operator=(const G4VEmProcess &right);

// =====================================================================

protected:

  G4ParticleChangeForGamma     fParticleChange;

private:

  G4EmModelManager*            modelManager;

  // tables and vectors
  G4PhysicsTable*              theLambdaTable;
  G4double*                    theEnergyOfCrossSectionMax;
  G4double*                    theCrossSectionMax;

  const G4ParticleDefinition*  particle;
  const G4ParticleDefinition*  secondaryParticle;
  const G4ParticleDefinition*  theGamma;
  const G4ParticleDefinition*  theElectron;
  const G4ParticleDefinition*  thePositron;

  const std::vector<G4double>* theCutsGamma;
  const std::vector<G4double>* theCutsElectron;
  const std::vector<G4double>* theCutsPositron;

  G4int                        nLambdaBins;

  G4double                     minKinEnergy;
  G4double                     maxKinEnergy;
  G4double                     lambdaFactor;

  // cash
  const G4Material*            currentMaterial;
  const G4MaterialCutsCouple*  currentCouple;
  size_t                       currentMaterialIndex;

  G4double                     mfpKinEnergy;
  G4double                     preStepKinEnergy;
  G4double                     preStepLambda;
  G4double                     preStepMFP;

  G4bool                       integral;
  G4bool                       meanFreePath;
  G4bool                       aboveCSmax;
  G4bool                       buildLambdaTable;
  G4bool                       applyCuts;
  G4bool                       startFromNull;

  G4int                        nRegions;
  std::vector<G4Region*>       regions;
  std::vector<G4bool>          flagsDeexcitation;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    if(!meanFreePath) ResetNumberOfInteractionLengthLeft();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambda(G4double& kineticEnergy,
                                  const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return GetCurrentLambda(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetCurrentLambda(G4double e)
{
  G4double x = 0.0;
  if(theLambdaTable) x = GetLambdaFromTable(e);
  else               x = ComputeCurrentLambda(e);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::RecalculateLambda(G4double e, const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return ComputeCurrentLambda(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::ComputeCurrentLambda(G4double e)
{
  G4VEmModel* currentModel = SelectModel(e);
  return currentModel->CrossSectionPerVolume(currentMaterial,particle,e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetLambdaFromTable(G4double e)
{
  G4bool b;
  return (((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::ComputeIntegralLambda(G4double e)
{
  meanFreePath  = false;
  aboveCSmax    = false;
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentMaterialIndex];
  if (e <= mfpKinEnergy) {
    preStepLambda = GetLambdaFromTable(e);
  } else {
    aboveCSmax  = true;
    G4double e1 = e*lambdaFactor;
    if(e1 > mfpKinEnergy) {
      preStepLambda  = GetLambdaFromTable(e);
      G4double preStepLambda1 = GetLambdaFromTable(e1);
      if(preStepLambda1 > preStepLambda) {
        mfpKinEnergy = e1;
        preStepLambda = preStepLambda1;
      }
    } else {
      preStepLambda = theCrossSectionMax[currentMaterialIndex];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetMeanFreePath(const G4Track& track,
                                                    G4double,
                                                    G4ForceCondition* condition)
{
  *condition = NotForced;
  preStepKinEnergy = track.GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  if( aboveCSmax && preStepKinEnergy < mfpKinEnergy ) ResetNumberOfInteractionLengthLeft();
  if (meanFreePath) {
    if(integral) ComputeIntegralLambda(preStepKinEnergy);
    else         preStepLambda = GetCurrentLambda(preStepKinEnergy);
    if(0.0 < preStepLambda) preStepMFP = 1.0/preStepLambda;
    else                    preStepMFP = DBL_MAX;
  }
  //    G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy<< " eCSmax= " 
  //  <<mfpKinEnergy<< " mfp= "<<preStepMFP<<G4endl;
  return preStepMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmProcess::SelectModel(G4double& kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmProcess::SelectModelForMaterial(
                                   G4double kinEnergy, size_t& idxRegion) const
{
  return modelManager->SelectModel(kinEnergy, idxRegion);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::ResetNumberOfInteractionLengthLeft()
{
  meanFreePath = true;
  aboveCSmax   = false;
  G4VProcess::ResetNumberOfInteractionLengthLeft();
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

inline G4double G4VEmProcess::GetGammaEnergyCut()
{
  return (*theCutsGamma)[currentMaterialIndex];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::GetElectronEnergyCut()
{
  return (*theCutsElectron)[currentMaterialIndex];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetLambdaFactor(G4double val)
{
  if(val > 0.0 && val <= 1.0) lambdaFactor = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
