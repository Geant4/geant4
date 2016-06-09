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
// $Id: G4VMscModel.hh,v 1.4 2008/03/10 10:39:28 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

private:

  //  hide assignment operator
  G4VMscModel & operator=(const  G4VMscModel &right);
  G4VMscModel(const  G4VMscModel&);

protected:

  G4double facrange;
  G4double facgeom;
  G4double facsafety;
  G4double skin;
  G4double dtrl;
  G4double lambdalimit;

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

#endif

