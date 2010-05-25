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
// ClassName:   G4MiscCHIPSBuilder
//
// Author: 2010, G.Folger
//
// created by modifying G4QInelasticCHIPSBuilder to apply only to MISC particles
// -----------------------------------------------------------------------------
#include "G4MiscCHIPSBuilder.hh"

G4MiscCHIPSBuilder::G4MiscCHIPSBuilder(G4int ver):
 verbose(ver), wasActivated(false)
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
}

G4MiscCHIPSBuilder::~G4MiscCHIPSBuilder()
{
  if(wasActivated) delete inelastic;
}

void G4MiscCHIPSBuilder::Build()
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
       pname == "pi-"   || pname == "pi+"   || pname == "neutron" ||  pname == "proton" )
    {
      if(verbose>1)G4cout<<"** G4MiscCHIPSBuilder: "<<pname<<" already defined"<<G4endl;
    }
    else  if( pname == "lambda"       || pname == "sigma+"       || pname == "sigma0"       ||
              pname == "sigma-"       || pname == "xi0" || pname == "xi-" || pname == "omega-" ||
              pname == "anti_proton"  || pname == "anti_neutron" || pname == "anti_lambda"  ||
              pname == "anti_sigma+"  || pname == "anti_sigma0"  || pname == "anti_sigma-"  ||
              pname == "anti_xi0"     || pname == "anti_xi-"     || pname == "anti_omega-"  )
    {
      if(verbose>1)G4cout<< "__G4QInelCHIPSBuilder: "<< pname <<" is defined here"<<G4endl;
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(inelastic);
      if(verbose>1) G4cout<<"###> G4MiscCHIPSBuilder: "<<inelastic->GetProcessName()
                          <<" is added for "<<pname<<G4endl;
    }
  }
}
