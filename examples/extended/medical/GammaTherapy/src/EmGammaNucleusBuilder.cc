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
// $Id: EmGammaNucleusBuilder.cc,v 1.6 2006/06/29 17:27:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#include "EmGammaNucleusBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4TheoFSGenerator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmGammaNucleusBuilder::EmGammaNucleusBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmGammaNucleusBuilder::~EmGammaNucleusBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmGammaNucleusBuilder::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Proton::Proton();
  G4Neutron::Neutron();
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmGammaNucleusBuilder::ConstructProcess()
{
  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4ProcessManager* pmanager = particle->GetProcessManager();
  G4PhotoNuclearProcess* pnp = new G4PhotoNuclearProcess("gNucler");
  G4TheoFSGenerator* tf = new G4TheoFSGenerator();
  G4GammaNuclearReaction* gn = new G4GammaNuclearReaction();
  tf->SetMinEnergy(3.5*GeV);
  tf->SetMaxEnergy(100.*TeV);
  pnp->RegisterMe(tf);
  gn->SetMaxEnergy(3.5*GeV);
  pnp->RegisterMe(gn);
  pmanager->AddDiscreteProcess(pnp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

