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
// $Id: G4QElasticPhysics.cc,v 1.2 2010-06-04 10:44:07 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QElasticPhysics
//
// Author: 3 Nov 2010 M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QElasticPhysics.hh"

#include "G4UHadronElasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronElastic.hh"
#include "G4QElastic.hh"

#include "G4VQCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

G4QElasticPhysics::G4QElasticPhysics(const G4String& name,  G4int ver)
  : G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  if(verbose > 1) G4cout << "### QElasticPhysics is initialized" << G4endl;
  model = 0;
}

void G4QElasticPhysics::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4QElasticPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  process = new G4QElastic();

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(pname == "anti_neutron" ||
       pname == "anti_proton"  ||
       pname == "anti_lambda"  ||
       pname == "anti_sigma-"  ||
       pname == "anti_sigma0"  ||
       pname == "anti_sigma+"  ||
       pname == "anti_xi-"     ||
       pname == "anti_xi0"     ||
       pname == "anti_omega-"  || 
       pname == "pi-"          ||
       pname == "pi+"          ||
       pname == "kaon-"        ||
       pname == "kaon+"        ||
       pname == "kaon0S"       ||
       pname == "kaon0L"       ||
       pname == "lambda"       ||
       pname == "sigma-"       ||
       pname == "sigma0"       ||
       pname == "sigma+"       ||
       pname == "xi-"          ||
       pname == "xi0"          ||
       pname == "omega-"       )
    {
      pmanager->AddDiscreteProcess(process);
      if(verbose>0) G4cout<<"### G4QElastic process is added for particle="<<pname<<G4endl;
    }
  }
}


