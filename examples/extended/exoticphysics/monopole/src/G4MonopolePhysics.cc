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
//
// $Id: G4MonopolePhysics.cc,v 1.1 2007/08/16 10:32:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4MonopolePhysics
//
// Author:      V.Ivanchenko 13.03.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Monopole.hh"
#include "G4MonopolePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4StepLimiter.hh"
#include "G4Transportation.hh"
#include "G4MultipleScattering.hh"
#include "G4mplIonisation.hh"
#include "G4hhIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysics::G4MonopolePhysics(const G4String& name):  G4VPhysicsConstructor(name)
{}

G4MonopolePhysics::~G4MonopolePhysics()
{}

void G4MonopolePhysics::ConstructParticle()
{
  G4Monopole::MonopoleDefinition(); 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4MonopolePhysics::ConstructProcess()
{

  // Add standard EM Processes for gamma
  G4cout << "G4MonopolePhysics::ConstructProcess" << G4endl;

  G4Monopole* mpl = G4Monopole::MonopoleDefinition();
  
  G4ProcessManager* pmanager = new G4ProcessManager(mpl);
  mpl->SetProcessManager(pmanager);
  
  pmanager->AddProcess( new G4Transportation(), -1, 0, 0);
  pmanager->AddProcess( new G4mplIonisation(mpl->MagneticCharge()), -1, 1, 1);
  pmanager->AddProcess( new G4StepLimiter(),  -1, -1, 3);
  if(mpl->GetPDGCharge() != 0.0) {
    pmanager->AddProcess(new G4hhIonisation(),  -1, 2, 2);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

