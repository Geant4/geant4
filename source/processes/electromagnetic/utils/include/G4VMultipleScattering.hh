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
// $Id: G4VMultipleScattering.hh 106714 2017-10-20 09:38:06Z gcosmo $
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
// 07-04-09 Moved msc methods from G4VEmModel to G4VMscModel (VI) 
//

// -------------------------------------------------------------------
//

#ifndef G4VMultipleScattering_h
#define G4VMultipleScattering_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4EmModelManager.hh"
#include "G4VMscModel.hh"
#include "G4EmParameters.hh"
#include "G4MscStepLimitType.hh"

class G4ParticleDefinition;
class G4VEnergyLossProcess;
class G4LossTableManager;
class G4SafetyHelper;

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

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override = 0;

  // obsolete
  virtual void PrintInfo() {};

  virtual void ProcessDescription(std::ostream& outFile) const override;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) = 0;

  virtual void StreamProcessInfo(std::ostream&, G4String) const {};

public:

  //------------------------------------------------------------------------
  // Generic methods common to all ContinuousDiscrete processes
  //------------------------------------------------------------------------

  // Initialise for build of tables
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // Build physics table during initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // Store PhysicsTable in a file.
  // Return false in case of failure at I/O
  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                           G4bool ascii = false) override;

  // Retrieve Physics from a file.
  // (return true if the Physics Table can be build by using file)
  // (return false if the process has no functionality or in case of failure)
  // File name should is constructed as processName+particleName and the
  // should be placed under the directory specifed by the argument.
  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                              G4bool ascii) override;

  // This is called in the beginning of tracking for a new track
  void StartTracking(G4Track*) override;

  // The function overloads the corresponding function of the base
  // class.It limits the step near to boundaries only
  // and invokes the method GetMscContinuousStepLimit at every step.
  G4double AlongStepGetPhysicalInteractionLength(
                                        const G4Track&,
                                        G4double  previousStepSize,
                                        G4double  currentMinimalStep,
                                        G4double& currentSafety,
                                        G4GPILSelection* selection) override;

  // The function overloads the corresponding function of the base
  // class.
  G4double PostStepGetPhysicalInteractionLength(
                                      const G4Track&,
                                      G4double  previousStepSize,
                                      G4ForceCondition* condition) override;

  // Along step actions
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  // Post step actions
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  // This method does not used for tracking, it is intended only for tests
  G4double ContinuousStepLimit(const G4Track& track,
                               G4double previousStepSize,
                               G4double currentMinimalStep,
                               G4double& currentSafety);

  //------------------------------------------------------------------------
  // Specific methods to set, access, modify models
  //------------------------------------------------------------------------

  // Select model in run time
  inline G4VEmModel* SelectModel(G4double kinEnergy, size_t idx);

public:

  // Add model for region, smaller value of order defines which
  // model will be selected for a given energy interval  
  void AddEmModel(G4int order, G4VEmModel*, const G4Region* region = nullptr);

  // Assign a model to a process local list, to enable the list in run time 
  // the derived process should execute AddEmModel(..) for all such models
  void SetEmModel(G4VMscModel*, size_t index = 0);
  
  // return a model from the local list
  G4VMscModel* EmModel(size_t index = 0) const;

  // Access to run time models by index
  inline G4VEmModel* GetModelByIndex(G4int idx = 0, G4bool ver = false) const;

  //------------------------------------------------------------------------
  // Get/Set parameters for simulation of multiple scattering
  //------------------------------------------------------------------------

  void SetIonisation(G4VEnergyLossProcess*);
  
  inline G4bool LateralDisplasmentFlag() const;
  inline void SetLateralDisplasmentFlag(G4bool val);

  inline G4double Skin() const;
  inline void SetSkin(G4double val);
  
  inline G4double RangeFactor() const;
  inline void SetRangeFactor(G4double val);
  
  inline G4double GeomFactor() const;
 
  inline G4double PolarAngleLimit() const;
 
  inline G4MscStepLimitType StepLimitType() const;
  inline void SetStepLimitType(G4MscStepLimitType val);

  inline G4double LowestKinEnergy() const;
  inline void SetLowestKinEnergy(G4double val);

  inline const G4ParticleDefinition* FirstParticle() const;

  //------------------------------------------------------------------------
  // Run time methods
  //------------------------------------------------------------------------

protected:

  // This method is not used for tracking, it returns mean free path value
  G4double GetMeanFreePath(const G4Track& track,
                           G4double,
                           G4ForceCondition* condition) override;

  // This method is not used for tracking, it returns step limit
  G4double GetContinuousStepLimit(const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimalStep,
                                  G4double& currentSafety) override ;

  // return number of already added models
  inline G4int NumberOfModels() const;

private:

  // hide  assignment operator
  G4VMultipleScattering(G4VMultipleScattering &) = delete;
  G4VMultipleScattering & 
    operator=(const G4VMultipleScattering &right) = delete;

  // Print out of generic class parameters
  void StreamInfo(std::ostream& outFile, const G4ParticleDefinition&,
                  G4String endOfLine=G4String("\n")) const;

  // ======== Parameters of the class fixed at construction =========

  G4EmModelManager*           modelManager;
  G4LossTableManager*         emManager;
  G4EmParameters*             theParameters;  

  // ======== Parameters of the class fixed at initialisation =======

  G4SafetyHelper*             safetyHelper;

  std::vector<G4VMscModel*>   mscModels;
  G4int                       numberOfModels;

  const G4ParticleDefinition* firstParticle;
  const G4ParticleDefinition* currParticle;

  G4MscStepLimitType          stepLimit;

  G4double                    facrange;
  G4double                    lowestKinEnergy;

  G4bool                      latDisplacement;
  G4bool                      isIon;
  G4bool                      fDispBeyondSafety;

  // ======== Cached values - may be state dependent ================

protected:

  G4ParticleChangeForMSC      fParticleChange;

private:

  // cache
  G4VMscModel*                currentModel;
  G4VEnergyLossProcess*       fIonisation;

  G4double                    geomMin;
  G4double                    minDisplacement2;
  G4double                    physStepLimit;
  G4double                    tPathLength;
  G4double                    gPathLength;

  G4ThreeVector               fNewPosition;
  G4ThreeVector               fNewDirection;
  G4bool                      fPositionChanged;
  G4bool                      isActive;
};

// ======== Run time inline methods ================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* 
G4VMultipleScattering::SelectModel(G4double kinEnergy, size_t coupleIndex)
{
  return modelManager->SelectModel(kinEnergy, coupleIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4bool G4VMultipleScattering::LateralDisplasmentFlag() const
{
  return latDisplacement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetLateralDisplasmentFlag(G4bool val)
{
  latDisplacement = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::Skin() const
{
  return theParameters->MscSkin();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetSkin(G4double val)
{
  theParameters->SetMscSkin(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::RangeFactor() const
{
  return facrange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  void G4VMultipleScattering::SetRangeFactor(G4double val)
{
  if(val > 0.0 && val < 1.0) { 
    facrange = val;
    theParameters->SetMscRangeFactor(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::GeomFactor() const
{
  return theParameters->MscGeomFactor();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline  G4double G4VMultipleScattering::PolarAngleLimit() const
{
  return theParameters->MscThetaLimit();
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
  if(val == fMinimal) { SetRangeFactor(0.2); }
  theParameters->SetMscStepLimitType(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering::LowestKinEnergy() const
{
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VMultipleScattering::SetLowestKinEnergy(G4double val)
{
  lowestKinEnergy = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4ParticleDefinition* G4VMultipleScattering::FirstParticle() const
{
  return firstParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VMultipleScattering::NumberOfModels() const
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* 
G4VMultipleScattering::GetModelByIndex(G4int idx, G4bool ver) const
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
