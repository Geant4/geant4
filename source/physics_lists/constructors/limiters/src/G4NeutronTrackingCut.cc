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
// $Id: G4NeutronTrackingCut.cc 102617 2017-02-10 07:58:33Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronTrackingCut
//
// Author: Nov 2006 G.Folger
//
//
//----------------------------------------------------------------------------
//

#include "G4NeutronTrackingCut.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Neutron.hh"
#include "G4NeutronKiller.hh"
#include "G4HadronicProcessStore.hh"

#include "G4Threading.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4NeutronTrackingCut);
//

G4NeutronTrackingCut::G4NeutronTrackingCut(G4int ver)
  :  G4VPhysicsConstructor("neutronTrackingCut")
   , verbose(ver)
{
  timeLimit          = 10.*microsecond;
  kineticEnergyLimit = 0.0;
}

G4NeutronTrackingCut::G4NeutronTrackingCut(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver)
{
  timeLimit          = 10.*microsecond;
  kineticEnergyLimit = 0.0;
}

G4NeutronTrackingCut::~G4NeutronTrackingCut()
{}

void G4NeutronTrackingCut::ConstructParticle()
{
  G4Neutron::NeutronDefinition();
}

void G4NeutronTrackingCut::ConstructProcess()
{
  G4NeutronKiller* pNeutronKiller = new G4NeutronKiller();
  G4ParticleDefinition * particle = G4Neutron::Neutron();
  G4ProcessManager * pmanager = particle->GetProcessManager();

  if(verbose > 0 && G4Threading::IsMasterThread()) {
    G4String pn = particle->GetParticleName();//Avoid data-race when passing
                                              //this string to G4MTcout
    G4cout << "### Adding tracking cuts for " << pn 
	   << "  TimeCut(ns)= " << timeLimit/ns 
	   << "  KinEnergyCut(MeV)= " <<  kineticEnergyLimit/MeV
	   <<  G4endl;
  }
  pmanager -> AddDiscreteProcess(pNeutronKiller);
  pNeutronKiller->SetKinEnergyLimit(kineticEnergyLimit);
  pNeutronKiller->SetTimeLimit(timeLimit);

  G4HadronicProcessStore::Instance()->RegisterExtraProcess(pNeutronKiller);
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(pNeutronKiller, particle);
}


