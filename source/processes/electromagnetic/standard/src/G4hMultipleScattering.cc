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
// $Id: G4hMultipleScattering.cc,v 1.3 2007/03/20 15:40:59 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-03 $
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
#include "G4UrbanMscModel.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4hMultipleScattering::G4hMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  lowKineticEnergy  = 0.1*keV;
  highKineticEnergy = 100.*TeV;
  totBins           = 120;

  facrange          = 0.2;
  dtrl              = 0.05;
  lambdalimit       = 1.*mm;
  facgeom           = 0.1;
  
  steppingAlgorithm = false;
  samplez           = false ; 
  isInitialized     = false;  

  SetBinning(totBins);
  SetMinKinEnergy(lowKineticEnergy);
  SetMaxKinEnergy(highKineticEnergy);

  SetLateralDisplasmentFlag(true);
  SetSkin(0.0);
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

void G4hMultipleScattering::MscStepLimitation(G4bool algorithm, G4double factor) 
{
  steppingAlgorithm = algorithm;
  if (factor > 0.) SetFacrange(factor);
  else { if (algorithm) SetFacrange(0.02); else SetFacrange(0.2);}

  if(verboseLevel > 1)  
    G4cout << "Stepping algorithm is set to " << steppingAlgorithm 
	   << " with facrange = " << facrange << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(isInitialized) {
    mscUrban->SetMscStepLimitation(steppingAlgorithm, facrange);
    if (p->GetParticleType() != "nucleus") {
      mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
      mscUrban->SetSkin(Skin());
    }
    return;
  }

  if (p->GetParticleType() == "nucleus") {
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
    SetSkin(0.0);
  } else {
    SetBuildLambdaTable(true);
  }
  mscUrban = new G4UrbanMscModel(facrange,dtrl,lambdalimit,
                                 facgeom,Skin(),
                                 samplez,steppingAlgorithm);
  mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  mscUrban->SetLowEnergyLimit(lowKineticEnergy);
  mscUrban->SetHighEnergyLimit(highKineticEnergy);
  AddEmModel(1,mscUrban);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::PrintInfo()
{
  G4cout << "      Boundary/stepping algorithm is active with facrange= "
	 << facrange
	 << "  Step limitation " << steppingAlgorithm
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

