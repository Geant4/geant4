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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QInelasticCHIPS_HPBuilder
//
// Author: 2009 M. Kosov
//
//----------------------------------------------------------------------------
//
// Short comment: This is a physics list of only one model G4QInelastic for
// all hadron-nuclear interactions at all energies + G4QNGamma for LE neutrons.
// There are process-mixings (G4QProcessMixer) w/ HP processes at low (< 20 MeV)
// energies with the G4QInelastic & G4QNGammafor neutrons, as CHIPS in G4 does
// not include yet all precise LE nA inelastic processes. In this particular
// builder the G4QInelastic process is attached to all hadrons + G4NGamma +
// HP_fission (fission && NNGamma CHIPS processes to be added).
//
// -----------------------------------------------------------------------------

//#define debug

#include "G4QInelasticCHIPS_HPBuilder.hh"
#include "G4SystemOfUnits.hh"

G4QInelasticCHIPS_HPBuilder::G4QInelasticCHIPS_HPBuilder(G4int ver):
    verbose(ver)
    , wasActivated(false)
    , inelastic(0)
    , theInProcessMixer(0)
    , theNgProcessMixer(0)
    , theFiProcessMixer(0) 
    , theNeutronInelastic(0)
    , theNeutronFission(0)
    , theNeutronCapture(0)
    , theCHIPSInelastic(0)
    , theCHIPSNGamma(0)
    , theHPNeutron(0)
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
}

G4QInelasticCHIPS_HPBuilder::~G4QInelasticCHIPS_HPBuilder()
{
  if(wasActivated)
  {
    delete inelastic;
    delete theCHIPSInelastic;
    delete theCHIPSNGamma;
    //delete theCHIPSFission;
    delete theNeutronInelastic;
    delete theNeutronCapture;
    delete theNeutronFission;
    delete theHPNeutron;
    delete theInProcessMixer;
    delete theNgProcessMixer;
    delete theFiProcessMixer;
  }
}

void G4QInelasticCHIPS_HPBuilder::Build()
{
  if(wasActivated) return;
  wasActivated = true;
  theParticleIterator->reset();
  inelastic = new G4QInelastic();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "kaon-" || pname == "kaon+" || pname == "kaon0S"  ||  pname == "kaon0L" ||
       //pname == "pi-" || pname == "pi+"   || pname == "neutron" ||  pname == "proton" ||
       pname == "pi-"   || pname == "pi+"   ||  pname == "proton" ||
       pname == "lambda"       || pname == "sigma+"       || pname == "sigma0"       ||
       pname == "sigma-"       || pname == "xi0" || pname == "xi-" || pname == "omega-" ||
       pname == "anti_proton"  || pname == "anti_neutron" || pname == "anti_lambda"  ||
       pname == "anti_sigma+"  || pname == "anti_sigma0"  || pname == "anti_sigma-"  ||
       pname == "anti_xi0"     || pname == "anti_xi-"     || pname == "anti_omega-"  )
    {
      if(verbose>1)
        G4cout<< "__G4QInelCHIPS_HPBuilder: "<< pname <<" is defined here"<<G4endl;
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(inelastic);
      if(verbose>1)G4cout<<"###>G4QInelasticCHIPS_HPBuilder: "<<inelastic->GetProcessName()
                         <<" is added for "<<pname<<G4endl;
    }
    else if(pname == "neutron")
    {
      if(verbose>1)
        G4cout<< "__G4QInelCHIPS_HPBuilder: "<< pname <<" is defined here"<<G4endl;
      G4ProcessManager* pmanager = particle->GetProcessManager();
      // The model definition for neutrons (needed for HP implementation)
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before NeutronBuild"<<G4endl;
#endif
      theCHIPSInelastic   = new G4QInelastic();
      theCHIPSNGamma      = new G4QNGamma();
      //theCHIPSFission   = new G4QFission();
      theNeutronInelastic = new G4NeutronInelasticProcess();
      theNeutronCapture   = new G4HadronCaptureProcess();
      theNeutronFission   = new G4HadronFissionProcess();
      theInProcessMixer   = new G4QDiscProcessMixer("Mixed NeutronInelastic", particle);
      theNgProcessMixer   = new G4QDiscProcessMixer("Mixed NGamma", particle);
      theFiProcessMixer   = new G4QDiscProcessMixer("Mixed NFission", particle);
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before Build HP processes"<<G4endl;
#endif
      theHPNeutron = new G4NeutronHPBuilder;
      theHPNeutron->Build(theNeutronInelastic);
      theHPNeutron->Build(theNeutronCapture);
      theHPNeutron->Build(theNeutronFission);
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before QIn="<<theCHIPSInelastic<<G4endl;
#endif
      theInProcessMixer->AddDiscreteProcess(theCHIPSInelastic, 1.E8*megaelectronvolt);
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: befr HPI="<<theNeutronInelastic<<G4endl;
#endif
      theInProcessMixer->AddDiscreteProcess(theNeutronInelastic, 19.9*megaelectronvolt);
      
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before QNG="<<theCHIPSNGamma<<G4endl;
#endif
      theNgProcessMixer->AddDiscreteProcess(theCHIPSNGamma, 1.E8*megaelectronvolt);
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before HPC="<<theNeutronCapture<<G4endl;
#endif
      theNgProcessMixer->AddDiscreteProcess(theNeutronCapture, 19.9*megaelectronvolt);
      
      //theFiProcessMixer->AddDiscreteProcess(theCHIPSFission, 1.E8*megaelectronvolt);
      //theFiProcessMixer->AddDiscreteProcess(theNeutronFission, 19.9*megaelectronvolt);
      
#ifdef debug
	G4cout<<"G4QInelasticCHIPS_HPBuilder::Build: before ProcessAdd"<<G4endl;
#endif
      pmanager->AddDiscreteProcess(theInProcessMixer); // Mix CHIPS+HP for neutronInelastic
      if(verbose>1)
        G4cout<<"###>G4QInelasticCHIPS_HPBuilder: "<<theCHIPSInelastic->GetProcessName()
              <<" is added for "<<pname<<G4endl;
      pmanager->AddDiscreteProcess(theNgProcessMixer);  // Mix CHIPS+HP for (n,gamma)
      if(verbose>1)
        G4cout<<"###>G4QInelasticCHIPS_HPBuilder: "<<theCHIPSNGamma->GetProcessName()
              <<" is added for "<<pname<<G4endl;
      pmanager->AddDiscreteProcess(theNeutronFission); // Only HP for fission
      //pmanager->AddDiscreteProcess(theFiProcessMixer); // Mix CHIPS+HP for fission
      if(verbose>1)
        G4cout<<"###>G4QInelasticCHIPS_HPBuilder: "<<theNeutronFission->GetProcessName()
              <<" is added for "<<pname<<G4endl;
    }
  }
}

// 2012 by M. Kosov
