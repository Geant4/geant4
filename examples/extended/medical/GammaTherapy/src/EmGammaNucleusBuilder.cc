//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// 
// $Id: EmGammaNucleusBuilder.cc,v 1.2 2004-12-03 11:06:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "EmGammaNucleusBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4GammaNuclearReaction.hh"

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
  G4GammaNuclearReaction* gn = new G4GammaNuclearReaction();
  gn->SetMaxEnergy(10.0*GeV);
  pnp->RegisterMe(gn);
  pmanager->AddDiscreteProcess(pnp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

