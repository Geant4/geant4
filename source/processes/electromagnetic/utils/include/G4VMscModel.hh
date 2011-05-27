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
// $Id: G4VMscModel.hh,v 1.10 2010-09-07 16:05:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//
// Class Description:
//
// General interface to msc models

// -------------------------------------------------------------------
//

#ifndef G4VMscModel_h
#define G4VMscModel_h 1

#include "G4VEmModel.hh"
#include "G4MscStepLimitType.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SafetyHelper.hh"
#include "G4VEnergyLossProcess.hh"

class G4ParticleChangeForMSC;

class G4VMscModel : public G4VEmModel
{

public:

  G4VMscModel(const G4String& nam);

  virtual ~G4VMscModel();

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,  
					      G4PhysicsTable*, 
					      G4double currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomPathLength);

  virtual void SampleScattering(const G4DynamicParticle*,
				G4double safety);

  // empty 
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double tmax);

  //================================================================
  //  Set parameters of multiple scattering models
  //================================================================
 
  inline void SetStepLimitType(G4MscStepLimitType);

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline void SetRangeFactor(G4double);

  inline void SetGeomFactor(G4double);

  inline void SetSkin(G4double);

  inline void SetSampleZ(G4bool);

  inline G4VEnergyLossProcess* GetIonisation() const;

  inline void SetIonisation(G4VEnergyLossProcess*);

protected:

  // initialisation of the ParticleChange for the model
  // initialisation of interface with geometry and ionisation 
  G4ParticleChangeForMSC* GetParticleChangeForMSC();

  // shift point of the track PostStep 
  void ComputeDisplacement(G4ParticleChangeForMSC*,  
			   const G4ThreeVector& displDir,
                           G4double displacement,
			   G4double postsafety);

  // compute safety
  inline G4double ComputeSafety(const G4ThreeVector& position, G4double limit);

  // compute linear distance to a geometry boundary
  inline G4double ComputeGeomLimit(const G4Track& position, G4double& presafety, 
				   G4double limit);

  inline G4double GetRange(const G4ParticleDefinition* part,
                           G4double kineticEnergy,
			   const G4MaterialCutsCouple* couple);

  inline G4double GetEnergy(const G4ParticleDefinition* part,
			    G4double range,
			    const G4MaterialCutsCouple* couple);

private:

  //  hide assignment operator
  G4VMscModel & operator=(const  G4VMscModel &right);
  G4VMscModel(const  G4VMscModel&);

  G4SafetyHelper* safetyHelper;
  G4VEnergyLossProcess* ionisation;

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

  G4MscStepLimitType steppingAlgorithm;

  G4bool   samplez;
  G4bool   latDisplasment;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetLateralDisplasmentFlag(G4bool val)
{
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetSkin(G4double val)
{
  skin = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetRangeFactor(G4double val)
{
  facrange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetGeomFactor(G4double val)
{
  facgeom = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetStepLimitType(G4MscStepLimitType val)
{
  steppingAlgorithm = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4VMscModel::SetSampleZ(G4bool val)
{
  samplez = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMscModel::ComputeSafety(const G4ThreeVector& position, 
					   G4double)
{
  return safetyHelper->ComputeSafety(position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4VMscModel::ComputeGeomLimit(const G4Track& track, 
					      G4double& presafety, 
					      G4double limit)
{
  G4double res = geomMax;
  if(track.GetVolume() != safetyHelper->GetWorldVolume()) {
    res = safetyHelper->CheckNextStep(
          track.GetStep()->GetPreStepPoint()->GetPosition(),
	  track.GetMomentumDirection(),
	  limit, presafety);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetRange(const G4ParticleDefinition* part,
		      G4double kinEnergy, const G4MaterialCutsCouple* couple)
{
  localrange = DBL_MAX;
  localtkin  = kinEnergy;
  if(ionisation) { localrange = ionisation->GetRangeForLoss(localtkin, couple); }
  else { 
    G4double q = part->GetPDGCharge()/eplus;
    localrange = kinEnergy/(dedx*q*q*couple->GetMaterial()->GetDensity()); 
  }
  return localrange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4VMscModel::GetEnergy(const G4ParticleDefinition* part,
		       G4double range, const G4MaterialCutsCouple* couple)
{
  G4double e = localtkin;
  if(ionisation) { e = ionisation->GetKineticEnergy(range, couple); }
  else { 
    G4double q = part->GetPDGCharge()/eplus;
    if(localrange > range) {
      e -= (localrange - range)*dedx*q*q*couple->GetMaterial()->GetDensity(); 
    }
  }
  if(e < 0.0) { e = 0.0; }
  return e;
}

inline G4VEnergyLossProcess* G4VMscModel::GetIonisation() const
{
  return ionisation;
}

inline void G4VMscModel::SetIonisation(G4VEnergyLossProcess* p)
{
  ionisation = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

