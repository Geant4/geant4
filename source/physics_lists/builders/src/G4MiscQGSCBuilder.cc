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
// ClassName:   G4MiscQGSCBuilder
//
// Author: 2009 M. Kosov (on the basis of the G4MiscLHEPBuilder)
//
//----------------------------------------------------------------------------
//
#include "G4MiscQGSCBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4MiscQGSCBuilder::G4MiscQGSCBuilder(G4int ver): 
    theModel(0),theCascade(0),theQGSCModel(0),
    theQGSCDecay(0),theQuasiElastic(0),
   verbose(ver), wasActivated(false)
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
}

G4MiscQGSCBuilder::~G4MiscQGSCBuilder() {}

void G4MiscQGSCBuilder::Build()
{
  wasActivated = true;

  //QGSC model definition
  theModel = new G4TheoFSGenerator("QGSC");  

  theQGSCModel    = new G4QGSModel< G4QGSParticipants >;
  theQGSCDecay    = new G4ExcitedStringDecay(new G4QGSMFragmentation);
  theQGSCModel->SetFragmentationModel(theQGSCDecay);
  theModel->SetHighEnergyGenerator(theQGSCModel);

  theCascade      = new G4StringChipsParticleLevelInterface;
  theModel->SetTransport(theCascade);

  theQuasiElastic = new G4QuasiElasticChannel;
  theModel->SetQuasiElasticChannel(theQuasiElastic);

  theModel->SetMinEnergy(0.);
  theModel->SetMaxEnergy(100*TeV);


  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "kaon-" || pname == "kaon+" || pname == "kaon0S"  ||  pname == "kaon0L" ||
       pname == "pi-"   || pname == "pi+"   || pname == "neutron" ||  pname == "proton" )
    {
      if(verbose>1)G4cout<<"** G4MiscQGSCBuilder: "<<pname<<" already defined"<<G4endl;
    }
    else if(
       pname == "anti_proton"  || pname == "anti_neutron" || pname == "anti_lambda"  ||
       pname == "anti_sigma+"  || pname == "anti_sigma0"  || pname == "anti_sigma-"  || 
       pname == "anti_xi0"     || pname == "anti_xi-"     || pname == "anti_omega-"  || 
       pname == "lambda"       || pname == "sigma+"       || pname == "sigma0"       ||
       pname == "sigma-"       || pname == "xi0" || pname == "xi-" || pname == "omega-")
     {
      if(verbose>1)G4cout<< "__ G4MiscQGSCBuilder: "<< pname <<" is defined here"<<G4endl;
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4HadronInelasticProcess* hp = new G4HadronInelasticProcess("hInelastic", particle);
      pmanager->AddDiscreteProcess(hp);
      hp->RegisterMe(theModel);
      if(verbose>1)
      G4cout<<"^^ G4MiscQGSCBuilder: "<<hp->GetProcessName()<<" added for "<<pname<<G4endl;
    }
  }
}

// 2009 by M. Kosov
