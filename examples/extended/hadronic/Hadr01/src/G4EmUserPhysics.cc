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
/// \file hadronic/Hadr01/src/G4EmUserPhysics.cc
/// \brief Implementation of the G4EmUserPhysics class
//
// $Id: G4EmUserPhysics.cc 70761 2013-06-05 12:30:51Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmUserPhysics
//
// Author:      V.Ivanchenko 11.07.2012
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4EmUserPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VProcess.hh"
#include "G4EmProcessOptions.hh"
#include "G4AntiProton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmUserPhysics::G4EmUserPhysics(G4int ver)
  : G4VPhysicsConstructor("User EM Options"), fVerbose(ver)
{
  G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmUserPhysics::~G4EmUserPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmUserPhysics::ConstructParticle()
{
  G4AntiProton::AntiProton();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmUserPhysics::ConstructProcess()
{
  const G4ParticleDefinition* pbar = G4AntiProton::AntiProton();
  G4ProcessManager* pmanager = pbar->GetProcessManager();
  G4ProcessVector* pv = pmanager->GetProcessList(); 
  size_t nn = pv->size();
  for(size_t i=0; i<nn; ++i) {
    G4VProcess* proc = (*pv)[i];
    if(fIonisation == proc->GetProcessSubType()) {
      G4VEnergyLossProcess* eloss = static_cast<G4VEnergyLossProcess*>(proc);
      G4double elim = 1.e-6*eV;
      eloss->SetLowestEnergyLimit(elim);
      G4cout << "### G4EmUserPhysics::ConstructProcess: "
             << "set new lowest energy limit " << elim/eV << " eV for "
             << pbar->GetParticleName() << G4endl;
      break;
    }
  }

  G4EmProcessOptions opt;
  opt.SetVerbose(fVerbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
