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
// $Id: G4VEmProcess.hh,v 1.47 2008/07/31 13:01:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  inline G4double RecalculateLambda(G4double kinEnergy,
				    const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Generic methods common to all Discrete processes 
  //------------------------------------------------------------------------

public:

  void PrintInfoDefinition();

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&);

  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&);

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

  //------------------------------------------------------------------------
  // Specific methods for Discrete EM post step simulation 
  //------------------------------------------------------------------------

  // It returns the cross section per volume for energy/ material
  G4double CrossSectionPerVolume(G4double kineticEnergy,
				 const G4MaterialCutsCouple* couple);

  // implementation of virtual method
  virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

  // It returns the cross section of the process per atom
  inline G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, 
					     G4double Z, G4double A=0., 
					     G4double cut=0.0);

  inline G4double MeanFreePath(const G4Track& track);

  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
					    size_t& idxRegion) const;

  // It returns cross section per volume
  inline G4double GetLambda(G4double& kinEnergy, 
			    const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Specific methods to build and access Physics Tables
  //------------------------------------------------------------------------

  // Binning for lambda table
  inline void SetLambdaBinning(G4int nbins);
  inline G4int LambdaBinning() const;

  // Min kinetic energy for tables
  inline void SetMinKinEnergy(G4double e);
  inline G4double MinKinEnergy() const;

  // Max kinetic energy for tables
  inline void SetMaxKinEnergy(G4double e);
  inline G4double MaxKinEnergy() const;

  inline void SetPolarAngleLimit(G4double a);
  inline G4double PolarAngleLimit() const;

  inline const G4PhysicsTable* LambdaTable() const;

  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

  inline const G4ParticleDefinition* Particle() const;
  inline const G4ParticleDefinition* SecondaryParticle() const;

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  // Add EM model coupled for the region
  inline void AddEmModel(G4int, G4VEmModel*, const G4Region* region = 0);
   
  // Assign a model to a process
  inline void SetModel(G4VEmModel*);
  
  // return the assigned model
  inline G4VEmModel* Model();
    
  // Define new energy range for the model identified by the name
  inline void UpdateEmModel(const G4String&, G4double, G4double);

  // Access to models
  inline G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false);

  //------------------------------------------------------------------------
  // Get/set parameters used for simulation of energy loss
  //------------------------------------------------------------------------

  inline void ActivateDeexcitation(G4bool, const G4Region* r = 0);

  inline void SetLambdaFactor(G4double val);

  inline void SetIntegral(G4bool val);
  inline G4bool IsIntegral() const;

  inline void SetApplyCuts(G4bool val);
  
protected:

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition);

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  inline G4ParticleChangeForGamma* GetParticleChange();

  inline void SetParticle(const G4ParticleDefinition* p);
  
  inline void SetSecondaryParticle(const G4ParticleDefinition* p);

  inline G4VEmModel* SelectModel(G4double& kinEnergy);

  inline size_t CurrentMaterialCutsCoupleIndex() const;

  inline G4double GetGammaEnergyCut();

  inline G4double GetElectronEnergyCut();

  inline void SetBuildTableFlag(G4bool val);

  inline void SetStartFromNullFlag(G4bool val);

private:

  void Clear();

  void BuildLambdaTable();

  void FindLambdaMax();

  inline void InitialiseStep(const G4Track&);

  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  inline void ComputeIntegralLambda(G4double kinEnergy);

  inline G4double GetLambdaFromTable(G4double kinEnergy);

  inline G4double GetCurrentLambda(G4double kinEnergy);

  inline G4double ComputeCurrentLambda(G4double kinEnergy);

  // hide  assignment operator

  G4VEmProcess(G4VEmProcess &);
  G4VEmProcess & operator=(const G4VEmProcess &right);

  // ======== Parameters of the class fixed at construction =========

  G4EmModelManager*            modelManager;
  const G4ParticleDefinition*  theGamma;
  const G4ParticleDefinition*  theElectron;
  const G4ParticleDefinition*  thePositron;
  const G4ParticleDefinition*  secondaryParticle;

  G4bool                       buildLambdaTable;

  // ======== Parameters of the class fixed at initialisation =======

  // tables and vectors
  G4PhysicsTable*              theLambdaTable;
  G4double*                    theEnergyOfCrossSectionMax;
  G4double*                    theCrossSectionMax;

  const std::vector<G4double>* theCuts;
  const std::vector<G4double>* theCutsGamma;
  const std::vector<G4double>* theCutsElectron;
  const std::vector<G4double>* theCutsPositron;

  G4int                        nLambdaBins;

  G4double                     minKinEnergy;
  G4double                     maxKinEnergy;
  G4double                     lambdaFactor;
  G4double                     polarAngleLimit;

  G4bool                       integral;
  G4bool                       applyCuts;
  G4bool                       startFromNull;

  G4int                        nRegions;
  std::vector<G4Region*>       regions;
  std::vector<G4bool>          flagsDeexcitation;

  // ======== Cashed values - may be state dependent ================

protected:

  G4ParticleChangeForGamma     fParticleChange;

private:

  std::vector<G4DynamicParticle*> secParticles;

  G4VEmModel*                  selectedModel;  

  const G4ParticleDefinition*  particle;

  // cash
  const G4Material*            currentMaterial;
  const G4MaterialCutsCouple*  currentCouple;
  size_t                       currentMaterialIndex;

  G4double                     mfpKinEnergy;
  G4double                     preStepKinEnergy;
  G4double                     preStepLambda;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    mfpKinEnergy = DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::InitialiseStep(const G4Track& track)
{
  preStepKinEnergy = track.GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  if (theNumberOfInteractionLengthLeft < 0.0) mfpKinEnergy = DBL_MAX;
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

inline G4double G4VEmProcess::RecalculateLambda(G4double e, 
						const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  return ComputeCurrentLambda(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::ComputeCurrentLambda(G4double e)
{
  G4VEmModel* currentModel = SelectModel(e);
  G4double x = 0.0;
  if(currentModel) 
    x = currentModel->CrossSectionPerVolume(currentMaterial,particle,
					    e,(*theCuts)[currentMaterialIndex]);
  return x;
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
  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentMaterialIndex];
  if (e <= mfpKinEnergy) {
    preStepLambda = GetLambdaFromTable(e);

  } else {
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

inline G4double G4VEmProcess::MeanFreePath(const G4Track& track)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepLambda = GetCurrentLambda(track.GetKineticEnergy());
  G4double x = DBL_MAX;
  if(DBL_MIN < preStepLambda) x = 1.0/preStepLambda;
  return x;
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

inline G4VEmModel* G4VEmProcess::GetModelByIndex(G4int idx, G4bool ver)
{
  return modelManager->GetModel(idx, ver);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::AddEmModel(G4int order, G4VEmModel* p, 
				     const G4Region* region)
{
  G4VEmFluctuationModel* fm = 0;
  modelManager->AddEmModel(order, p, fm, region);
  if(p) p->SetParticleChange(pParticleChange);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetModel(G4VEmModel* model)
{
  selectedModel = model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmProcess::Model()
{
  return selectedModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::UpdateEmModel(const G4String& nam, 
					G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::ComputeCrossSectionPerAtom(
		 G4double kineticEnergy, G4double Z, G4double A, G4double cut)
{
  G4VEmModel* model = SelectModel(kineticEnergy);
  return model->ComputeCrossSectionPerAtom(particle,kineticEnergy,Z,A,cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetLambdaBinning(G4int nbins)
{
  nLambdaBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmProcess::LambdaBinning() const
{
  return nLambdaBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetPolarAngleLimit(G4double val)
{
  if(val < 0.0)     polarAngleLimit = 0.0;
  else if(val > pi) polarAngleLimit = pi;
  else              polarAngleLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmProcess::PolarAngleLimit() const
{
  return polarAngleLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::ActivateDeexcitation(G4bool, const G4Region*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4PhysicsTable* G4VEmProcess::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetIntegral(G4bool val)
{
  if(particle && particle != theGamma) integral = val;
  if(integral) buildLambdaTable = true;
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
  if(!val) integral = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetStartFromNullFlag(G4bool val)
{
  startFromNull = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmProcess::SetApplyCuts(G4bool val)
{
  applyCuts = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline size_t G4VEmProcess::CurrentMaterialCutsCoupleIndex() const 
{
  return currentMaterialIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
