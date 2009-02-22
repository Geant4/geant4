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
// $Id: G4VMscModel.hh,v 1.7 2009-02-22 17:32:08 vnivanch Exp $
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
//
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

class G4ParticleChangeForMSC;

class G4VMscModel : public G4VEmModel
{

public:

  G4VMscModel(const G4String& nam);

  virtual ~G4VMscModel();

  inline void SetStepLimitType(G4MscStepLimitType);

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline void SetRangeFactor(G4double);

  inline void SetGeomFactor(G4double);

  inline void SetSkin(G4double);

  inline void SetSampleZ(G4bool);

protected:

  void InitialiseSafetyHelper();

  void ComputeDisplacement(G4ParticleChangeForMSC*,  
			   const G4ThreeVector& displDir,
                           G4double displacement,
			   G4double postsafety);

  inline G4double ComputeSafety(const G4ThreeVector& position, G4double limit);

  inline G4double ComputeGeomLimit(const G4Track& position, G4double& presafety, 
				   G4double limit);

private:

  //  hide assignment operator
  G4VMscModel & operator=(const  G4VMscModel &right);
  G4VMscModel(const  G4VMscModel&);

  G4SafetyHelper* safetyHelper;

protected:

  G4double facrange;
  G4double facgeom;
  G4double facsafety;
  G4double skin;
  G4double dtrl;
  G4double lambdalimit;
  G4double geommax;

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
  G4double res = geommax;
  if(track.GetVolume() != safetyHelper->GetWorldVolume()) {
    res = safetyHelper->CheckNextStep(
          track.GetStep()->GetPreStepPoint()->GetPosition(),
	  track.GetMomentumDirection(),
	  limit, presafety);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

