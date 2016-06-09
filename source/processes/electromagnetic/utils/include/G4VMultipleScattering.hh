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
// $Id: G4VMultipleScattering.hh,v 1.54 2008/07/31 13:01:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VMultipleScattering
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 12.03.2002
//
// Modifications:
//
// 16-07-03 Update GetRange interface (V.Ivanchenko)
//
//
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles
//
// 26-11-03 bugfix in AlongStepDoIt (L.Urban)
// 25-05-04 add protection against case when range is less than steplimit (VI)
// 30-06-04 make destructor virtual (V.Ivanchenko)
// 27-08-04 Add InitialiseForRun method (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 15-04-05 optimize internal interfaces (V.Ivanchenko)
// 15-04-05 remove boundary flag (V.Ivanchenko)
// 07-10-05 error in a protection in GetContinuousStepLimit corrected (L.Urban)
// 27-10-05 introduce virtual function MscStepLimitation() (V.Ivanchenko)
// 26-01-06 Rename GetRange -> GetRangeFromRestricteDEDX (V.Ivanchenko)
// 17-02-06 Save table of transport cross sections not mfp (V.Ivanchenko)
// 07-03-06 Move step limit calculation to model (V.Ivanchenko)
// 13-05-06 Add method to access model by index (V.Ivanchenko)
// 12-02-07 Add get/set skin (V.Ivanchenko)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
// 15-07-08 Reorder class members for further multi-thread development (VI)
//

// -------------------------------------------------------------------
//

#ifndef G4VMultipleScattering_h
#define G4VMultipleScattering_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4EmModelManager.hh"
#include "G4VEmModel.hh"
#include "G4MscStepLimitType.hh"

class G4ParticleDefinition;
class G4DataVector;
class G4PhysicsTable;
class G4PhysicsVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VMultipleScattering : public G4VContinuousDiscreteProcess
{
public:

  G4VMultipleScattering(const G4String& name = "msc",
			G4ProcessType type = fElectromagnetic);

  virtual ~G4VMultipleScattering();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) = 0;

  virtual void PrintInfo() = 0;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

public:

  //------------------------------------------------------------------------
  // Generic methods common to all ContinuousDiscrete processes
  //------------------------------------------------------------------------

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&);
  
  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&);

  // Print out of generic class parameters
  void PrintInfoDefinition();

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

  //------------------------------------------------------------------------
  // Specific methods for msc processes
  //------------------------------------------------------------------------

  // The function overloads the corresponding function of the base
  // class.It limits the step near to boundaries only
  // and invokes the method GetMscContinuousStepLimit at every step.
  virtual G4double AlongStepGetPhysicalInteractionLength(
                                            const G4Track&,
					    G4double  previousStepSize,
					    G4double  currentMinimalStep,
					    G4double& currentSafety,
					    G4GPILSelection* selection);

  // The function overloads the corresponding function of the base
  // class.
  G4double PostStepGetPhysicalInteractionLength(
                                            const G4Track&,
					    G4double  previousStepSize,
					    G4ForceCondition* condition);

  // This method does not used for tracking, it is intended only for tests
  inline G4double ContinuousStepLimit(const G4Track& track,
				      G4double previousStepSize,
				      G4double currentMinimalStep,
				      G4double& currentSafety);

  //------------------------------------------------------------------------
  // Specific methods to build and access Physics Tables
  //------------------------------------------------------------------------

  // Build empty Physics Vector
  G4PhysicsVector* PhysicsVector(const G4MaterialCutsCouple*);

  inline void SetBinning(G4int nbins);
  inline G4int Binning() const;

  inline void SetMinKinEnergy(G4double e);
  inline G4double MinKinEnergy() const;

  inline void SetMaxKinEnergy(G4double e);
  inline G4double MaxKinEnergy() const;

  inline void SetBuildLambdaTable(G4bool val);

  inline G4PhysicsTable* LambdaTable() const;

  //------------------------------------------------------------------------
  // Define and access particle type 
  //------------------------------------------------------------------------

  inline const G4ParticleDefinition* Particle() const;
  inline void SetParticle(const G4ParticleDefinition*);

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  inline void AddEmModel(G4int, G4VEmModel*, const G4Region* region = 0);

  inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy, 
					    size_t& idxRegion) const;

  // Access to models
  inline G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false);

  //------------------------------------------------------------------------
  // Set parameters for simulation of multiple scattering
  //------------------------------------------------------------------------

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline void SetSkin(G4double val);

  inline void SetRangeFactor(G4double val);

  inline void SetGeomFactor(G4double val);

  inline void SetPolarAngleLimit(G4double val);

  inline void SetStepLimitType(G4MscStepLimitType val);

protected:

  // This method is used for tracking, it returns mean free path value
  G4double GetMeanFreePath(const G4Track& track,
			   G4double,
			   G4ForceCondition* condition);

  //------------------------------------------------------------------------
  // Run time methods
  //------------------------------------------------------------------------

  // This method is not used for tracking, it returns step limit
  G4double GetContinuousStepLimit(const G4Track& track,
				  G4double previousStepSize,
				  G4double currentMinimalStep,
				  G4double& currentSafety);

  inline G4double GetLambda(const G4ParticleDefinition* p, 
			    G4double& kineticEnergy);

  // This method is used for tracking, it returns step limit
  inline G4double GetMscContinuousStepLimit(const G4Track& track,
					    G4double scaledKinEnergy,
					    G4double currentMinimalStep,
					    G4double& currentSafety);

  inline G4VEmModel* SelectModel(G4double kinEnergy);
  // Select concrete model

  inline const G4MaterialCutsCouple* CurrentMaterialCutsCouple() const; 

  // define current material
  inline void DefineMaterial(const G4MaterialCutsCouple* couple);

  //------------------------------------------------------------------------
  // Access parameters of multiple scattering
  //------------------------------------------------------------------------

  inline G4ParticleChangeForMSC* GetParticleChange();

  inline G4double Skin() const;

  inline G4double RangeFactor() const;

  inline G4double GeomFactor() const;

  inline G4double PolarAngleLimit() const;

  inline G4MscStepLimitType StepLimitType() const;

  inline G4bool LateralDisplasmentFlag() const;

private:

  // hide  assignment operator

  G4VMultipleScattering(G4VMultipleScattering &);
  G4VMultipleScattering & operator=(const G4VMultipleScattering &right);

  // ======== Parameters of the class fixed at construction =========

  G4EmModelManager*           modelManager;
  G4bool                      buildLambdaTable;

  // ======== Parameters of the class fixed at initialisation =======

  G4PhysicsTable*             theLambdaTable;
  const G4ParticleDefinition* firstParticle;

  G4MscStepLimitType          stepLimit;

  G4double                    minKinEnergy;
  G4double                    maxKinEnergy;
  G4double                    skin;
  G4double                    facrange;
  G4double                    facgeom;
  G4double                    polarAngleLimit;

  G4int                       nBins;

  G4bool                      latDisplasment;

  // ======== Cashed values - may be state dependent ================

protected:

  G4GPILSelection             valueGPILSelectionMSC;
  G4ParticleChangeForMSC      fParticleChange;

private:

  G4VEmModel*                 currentModel;

  // cache
  const G4ParticleDefinition* currentParticle;
  const G4MaterialCutsCouple* currentCouple;
  size_t                      currentMaterialIndex;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4VMultipleScattering::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterialIndex = couple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::GetMscContinuousStepLimit(
                                          const G4Track& track,
					  G4double scaledKinEnergy,
					  G4double currentMinimalStep,
					  G4double&)
{
  G4double x = currentMinimalStep;
  DefineMaterial(track.GetMaterialCutsCouple());
  currentModel = SelectModel(scaledKinEnergy);
  if(x > 0.0 && scaledKinEnergy > 0.0) {
    G4double tPathLength = 
      currentModel->ComputeTruePathLengthLimit(track, theLambdaTable, x);
    if (tPathLength < x) valueGPILSelectionMSC = CandidateForSelection;
    x = currentModel->ComputeGeomPathLength(tPathLength);  
    //  G4cout << "tPathLength= " << tPathLength
    //         << " stepLimit= " << x 
    //        << " currentMinimalStep= " << currentMinimalStep<< G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::ContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  return GetMscContinuousStepLimit(track,previousStepSize,currentMinimalStep,
				   currentSafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4VMultipleScattering::GetLambda(const G4ParticleDefinition* p, 
					  G4double& e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = currentModel->CrossSection(currentCouple,p,e);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX; 
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VMultipleScattering::SelectModel(G4double kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VMultipleScattering::SelectModelForMaterial(
		   G4double kinEnergy, size_t& idxRegion) const
{
  return modelManager->SelectModel(kinEnergy, idxRegion);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetBinning(G4int nbins)
{
  nBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VMultipleScattering::Binning() const
{
  return nBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4bool G4VMultipleScattering::LateralDisplasmentFlag() const
{
  return latDisplasment;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetLateralDisplasmentFlag(G4bool val)
{
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4ParticleChangeForMSC* G4VMultipleScattering::GetParticleChange()
{
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::Skin() const
{
  return skin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetSkin(G4double val)
{
  if(val < 1.0) skin = 0.0;
  else          skin = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::RangeFactor() const
{
  return facrange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetRangeFactor(G4double val)
{
  if(val > 0.0) facrange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::GeomFactor() const
{
  return facgeom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetGeomFactor(G4double val)
{
  if(val > 0.0) facgeom = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::PolarAngleLimit() const
{
  return polarAngleLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetPolarAngleLimit(G4double val)
{
  if(val < 0.0)     polarAngleLimit = 0.0;
  else if(val > pi) polarAngleLimit = pi;
  else              polarAngleLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4MscStepLimitType G4VMultipleScattering::StepLimitType() const
{
  return stepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetStepLimitType(G4MscStepLimitType val) 
{
  stepLimit = val;
  if(val == fMinimal) facrange = 0.2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetBuildLambdaTable(G4bool val)
{
  buildLambdaTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  const G4ParticleDefinition* G4VMultipleScattering::Particle() const
{
  return currentParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VMultipleScattering::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
const G4MaterialCutsCouple* G4VMultipleScattering::CurrentMaterialCutsCouple() const
{
  return currentCouple;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::AddEmModel(G4int order, G4VEmModel* p,
					      const G4Region* region)
{
  G4VEmFluctuationModel* fm = 0;
  modelManager->AddEmModel(order, p, fm, region);
  if(p) p->SetParticleChange(pParticleChange);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4VMultipleScattering::GetModelByIndex(G4int idx, G4bool ver)
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
