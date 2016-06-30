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
// $Id: G4VMscModel.hh 96626 2016-04-27 08:36:27Z gcosmo $
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
#include "G4VEnergyLossProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4ParticleChangeForMSC;

class G4VMscModel : public G4VEmModel
{

public:

  explicit G4VMscModel(const G4String& nam);

  virtual ~G4VMscModel();

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,  
					      G4double& stepLimit);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomPathLength);

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&,
					  G4double safety);

  // empty method
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double tmax) override;

  //================================================================
  //  Set parameters of multiple scattering models
  //================================================================
 
  inline void SetStepLimitType(G4MscStepLimitType);

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline void SetRangeFactor(G4double);

  inline void SetGeomFactor(G4double);

  inline void SetSkin(G4double);

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

public:

  // compute safety
  inline G4double ComputeSafety(const G4ThreeVector& position, 
				G4double limit= DBL_MAX);

  // compute linear distance to a geometry boundary
  inline G4double ComputeGeomLimit(const G4Track&, G4double& presafety, 
				   G4double limit);

  inline G4double GetDEDX(const G4ParticleDefinition* part,
			  G4double kineticEnergy,
			  const G4MaterialCutsCouple* couple);

  inline G4double GetRange(const G4ParticleDefinition* part,
                           G4double kineticEnergy,
			   const G4MaterialCutsCouple* couple);

  inline G4double GetEnergy(const G4ParticleDefinition* part,
			    G4double range,
			    const G4MaterialCutsCouple* couple);

  // G4MaterialCutsCouple should be defined before call to this method
  inline 
  G4double GetTransportMeanFreePath(const G4ParticleDefinition* part,
				    G4double kinEnergy);

private:

  //  hide assignment operator
  G4VMscModel & operator=(const  G4VMscModel &right) = delete;
  G4VMscModel(const  G4VMscModel&) = delete;

  G4SafetyHelper* safetyHelper;
  G4VEnergyLossProcess* ionisation;
  const G4ParticleDefinition* currentPart;

  G4double dedx;
  G4double localtkin;
  G4double localrange;

protected:

  G4double facrange;
  G4double facgeom;
  G4double facsafety;
  G4double skin;
  G4double dtrl;
  G4double lambdalimit;
  G4double geomMin;
  G4double geomMax;

  G4ThreeVector      fDisplacement;
  G4MscStepLimitType steppingAlgorithm;

  G4bool   samplez;
  G4bool   latDisplasment;

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

inline G4double 
G4VMscModel::GetDEDX(const G4ParticleDefinition* part,
		     G4double kinEnergy, const G4MaterialCutsCouple* couple)
{
  G4double x;
  if(ionisation) { x = ionisation->GetDEDX(kinEnergy, couple); }
  else { 
    G4double q = part->GetPDGCharge()*inveplus;
    x = dedx*q*q;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetRange(const G4ParticleDefinition* part,
		      G4double kinEnergy, const G4MaterialCutsCouple* couple)
{
  //G4cout << "G4VMscModel::GetRange E(MeV)= " << kinEnergy << "  " 
  //  << ionisation << "  " << part->GetParticleName()
  //	 << G4endl;
  localtkin  = kinEnergy;
  if(ionisation) { 
    localrange = ionisation->GetRangeForLoss(kinEnergy, couple); 
  } else { 
    G4double q = part->GetPDGCharge()*inveplus;
    localrange = kinEnergy/(dedx*q*q*couple->GetMaterial()->GetDensity()); 
  }
  //G4cout << "R(mm)= " << localrange << "  "  << ionisation << G4endl;
  return localrange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetEnergy(const G4ParticleDefinition* part,
		       G4double range, const G4MaterialCutsCouple* couple)
{
  G4double e;
  //G4cout << "G4VMscModel::GetEnergy R(mm)= " << range << "  " << ionisation
  //	 << "  Rlocal(mm)= " << localrange << "  Elocal(MeV)= " << localtkin
  //	 << G4endl;
  if(ionisation) { e = ionisation->GetKineticEnergy(range, couple); }
  else { 
    e = localtkin;
    if(localrange > range) {
      G4double q = part->GetPDGCharge()*inveplus;
      e -= (localrange - range)*dedx*q*q*couple->GetMaterial()->GetDensity(); 
    } 
  }
  return e;
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
  if(xSectionTable) {
    G4int idx = CurrentCouple()->GetIndex();
    x = (*xSectionTable)[(*theDensityIdx)[idx]]->Value(ekin, idxTable)
      *(*theDensityFactor)[idx]/(ekin*ekin);
  } else { 
    x = CrossSectionPerVolume(CurrentCouple()->GetMaterial(), part, ekin, 
			      0.0, DBL_MAX); 
  }
  x = (x > 0.0) ? 1.0/x : DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

