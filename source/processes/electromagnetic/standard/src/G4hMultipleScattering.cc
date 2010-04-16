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
// $Id: G4hMultipleScattering.cc,v 1.18 2010-04-16 09:01:24 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4hMultipleScattering
//
// Author:        Laszlo Urban
//
// Creation date: 24.10.2006 cloned from G4MultipleScattering
// 
// Modified:
// 12-02-07 skin can be changed via UI command (VI)
// 20.03.07 Remove local parameter skin, set facgeom=0.1(V.Ivanchenko) 
//
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4hMultipleScattering.hh"
#include "G4UrbanMscModel90.hh"
#include "G4MscStepLimitType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4hMultipleScattering::G4hMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  isInitialized = false;  
  SetStepLimitType(fMinimal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4hMultipleScattering::~G4hMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4hMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  if(isInitialized) return;
  G4VMscModel* mscUrban = new G4UrbanMscModel90();
  AddEmModel(1,mscUrban);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::PrintInfo()
{
  G4cout << "      RangeFactor= " << RangeFactor()
	 << ", stepLimitType: " << StepLimitType()
         << ", latDisplacement: " << LateralDisplasmentFlag();
  if(StepLimitType() == fUseDistanceToBoundary) {
    G4cout  << ", skin= " << Skin() << ", geomFactor= " << GeomFactor();
  }  
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

