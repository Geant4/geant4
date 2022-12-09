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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VMscModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2008
//
// Modifications:
// 07.04.2009 V.Ivanchenko moved msc methods from G4VEmModel to G4VMscModel 
// 26.03.2012 V.Ivanchenko added transport x-section pointer and Get?Set methods
//
// Class Description:
//
// General interface to msc models

// -------------------------------------------------------------------
//
#ifndef G4VMscModel_h
#define G4VMscModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4MscStepLimitType.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SafetyHelper.hh"
#include "G4PhysicsTable.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4ParticleChangeForMSC;
class G4ParticleDefinition;
class G4VEnergyLossProcess;

class G4VMscModel : public G4VEmModel
{

public:

  explicit G4VMscModel(const G4String& nam);

  ~G4VMscModel() override;

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,  
					      G4double& stepLimit) = 0;

  virtual G4double ComputeGeomPathLength(G4double truePathLength) = 0;

  virtual G4double ComputeTrueStepLength(G4double geomPathLength) = 0;

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&,
					  G4double safety) = 0;

  void InitialiseParameters(const G4ParticleDefinition*);

  void DumpParameters(std::ostream& out) const;

  // empty method
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin, G4double tmax) override;

  //================================================================
  //  Set parameters of multiple scattering models
  //================================================================
 
  inline void SetStepLimitType(G4MscStepLimitType);

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline void SetRangeFactor(G4double);

  inline void SetGeomFactor(G4double);

  inline void SetSkin(G4double);

  inline void SetLambdaLimit(G4double);

  inline void SetSafetyFactor(G4double);

  inline void SetSampleZ(G4bool);

  //================================================================
  //  Get/Set access to Physics Tables
  //================================================================

  inline G4VEnergyLossProcess* GetIonisation() const;

  inline void SetIonisation(G4VEnergyLossProcess*, 
			    const G4ParticleDefinition* part);

  //================================================================
  //  Run time methods
  //================================================================

protected:

  // initialisation of the ParticleChange for the model
  // initialisation of interface with geometry and ionisation 
  G4ParticleChangeForMSC* 
  GetParticleChangeForMSC(const G4ParticleDefinition* p = nullptr);

  // convert true length to geometry length
  inline G4double ConvertTrueToGeom(G4double& tLength, G4double& gLength);

  // should be set before initialisation
  inline void SetUseSplineForMSC(G4bool val);

public:

  // compute safety
  inline G4double ComputeSafety(const G4ThreeVector& position, 
				G4double limit= DBL_MAX);

  // compute linear distance to a geometry boundary
  inline G4double ComputeGeomLimit(const G4Track&, G4double& presafety, 
				   G4double limit);

  G4double GetDEDX(const G4ParticleDefinition* part,
                          G4double kineticEnergy,
                          const G4MaterialCutsCouple* couple);

  G4double GetDEDX(const G4ParticleDefinition* part,
                          G4double kineticEnergy,
                          const G4MaterialCutsCouple* couple,
                          G4double logKineticEnergy);

  G4double GetRange(const G4ParticleDefinition* part,
                           G4double kineticEnergy,
                           const G4MaterialCutsCouple* couple);

  G4double GetRange(const G4ParticleDefinition* part,
                           G4double kineticEnergy,
                           const G4MaterialCutsCouple* couple,
                           G4double logKineticEnergy);

  G4double GetEnergy(const G4ParticleDefinition* part,
			    G4double range,
			    const G4MaterialCutsCouple* couple);

  // G4MaterialCutsCouple should be defined before call to this method
  inline 
  G4double GetTransportMeanFreePath(const G4ParticleDefinition* part,
                                    G4double kinEnergy);

  inline
  G4double GetTransportMeanFreePath(const G4ParticleDefinition* part,
                                    G4double kinEnergy,
                                    G4double logKinEnergy);

  //  hide assignment operator
  G4VMscModel & operator=(const  G4VMscModel &right) = delete;
  G4VMscModel(const  G4VMscModel&) = delete;

private:

  G4SafetyHelper* safetyHelper = nullptr;
  G4VEnergyLossProcess* ionisation = nullptr;
  const G4ParticleDefinition* currentPart = nullptr;

  G4double dedx = 0.0;
  G4double localtkin = 0.0;
  G4double localrange = DBL_MAX;

protected:

  G4double facrange = 0.04;
  G4double facgeom = 2.5;
  G4double facsafety = 0.6;
  G4double skin = 1.0;
  G4double dtrl = 0.05;
  G4double lambdalimit;
  G4double geomMin;
  G4double geomMax;

  G4ThreeVector fDisplacement;
  G4MscStepLimitType steppingAlgorithm;

  G4bool samplez = false;
  G4bool latDisplasment = true;

private:

  G4bool useSpline = true;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetLateralDisplasmentFlag(G4bool val)
{
  if(!IsLocked()) { latDisplasment = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetSkin(G4double val)
{
  if(!IsLocked()) { skin = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetRangeFactor(G4double val)
{
  if(!IsLocked()) { facrange = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetGeomFactor(G4double val)
{
  if(!IsLocked()) { facgeom = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetLambdaLimit(G4double val)
{
  if(!IsLocked()) { lambdalimit = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetSafetyFactor(G4double val)
{
  if(!IsLocked()) { facsafety = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetStepLimitType(G4MscStepLimitType val)
{
  if(!IsLocked()) { steppingAlgorithm = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetSampleZ(G4bool val)
{
  if(!IsLocked()) { samplez = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMscModel::ComputeSafety(const G4ThreeVector& position, 
					   G4double limit)
{
   return safetyHelper->ComputeSafety(position, limit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMscModel::ConvertTrueToGeom(G4double& tlength, 
					       G4double& glength)
{
  glength = ComputeGeomPathLength(tlength);
  // should return true length 
  return tlength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMscModel::ComputeGeomLimit(const G4Track& track, 
					      G4double& presafety, 
					      G4double limit)
{
  return safetyHelper->CheckNextStep(
          track.GetStep()->GetPreStepPoint()->GetPosition(),
	  track.GetMomentumDirection(),
	  limit, presafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4VEnergyLossProcess* G4VMscModel::GetIonisation() const
{
  return ionisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetIonisation(G4VEnergyLossProcess* p,
				       const G4ParticleDefinition* part)
{
  ionisation = p;
  currentPart = part;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetTransportMeanFreePath(const G4ParticleDefinition* part,
                                      G4double ekin)
{
  G4double x;
  if (nullptr != xSectionTable) {
    x = pFactor*(*xSectionTable)[basedCoupleIndex]->Value(ekin)/(ekin*ekin);
  } else { 
    x = pFactor*CrossSectionPerVolume(pBaseMaterial, part, ekin, 0.0, DBL_MAX); 
  }
  return (x > 0.0) ? 1.0/x : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetTransportMeanFreePath(const G4ParticleDefinition* part,
                                      G4double ekin, G4double logekin)
{
  G4double x;
  if (nullptr != xSectionTable) {
    x = pFactor*(*xSectionTable)[basedCoupleIndex]->LogVectorValue(ekin, logekin)/(ekin*ekin);
  } else { 
    x = pFactor*CrossSectionPerVolume(pBaseMaterial, part, ekin, 0.0, DBL_MAX);
  }
  return (x > 0.0) ? 1.0/x : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetUseSplineForMSC(G4bool val)
{
  useSpline = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
