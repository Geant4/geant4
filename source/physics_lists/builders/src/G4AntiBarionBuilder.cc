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
// $Id: G4AntiBarionBuilder.cc 103801 2017-04-27 13:59:03Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4AntiBarionBuilder
//
// Author: 2011 J. Apostolakis
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#include "G4AntiBarionBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4AntiBarionBuilder::
G4AntiBarionBuilder(): wasActivated(false) 
{  
  theAntiProtonInelastic=new G4AntiProtonInelasticProcess;
  theAntiNeutronInelastic=new G4AntiNeutronInelasticProcess;
  theAntiDeuteronInelastic=new G4AntiDeuteronInelasticProcess;
  theAntiTritonInelastic=new G4AntiTritonInelasticProcess;
  theAntiHe3Inelastic=new G4AntiHe3InelasticProcess;
  theAntiAlphaInelastic=new G4AntiAlphaInelasticProcess;
}


void G4AntiBarionBuilder::
Build()
{
  wasActivated = true;

  std::vector<G4VAntiBarionBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theAntiProtonInelastic);
    (*i)->Build(theAntiNeutronInelastic);
    (*i)->Build(theAntiDeuteronInelastic);
    (*i)->Build(theAntiTritonInelastic);
    (*i)->Build(theAntiHe3Inelastic);
    (*i)->Build(theAntiAlphaInelastic);
  }
  G4ProcessManager * theProcMan;
  theProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiProtonInelastic);
  
  theProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiNeutronInelastic);
  
  theProcMan = G4AntiDeuteron::AntiDeuteron()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiDeuteronInelastic);

  theProcMan = G4AntiTriton::AntiTriton()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiTritonInelastic);
  
  theProcMan = G4AntiHe3::AntiHe3()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiHe3Inelastic);
  
  theProcMan = G4AntiAlpha::AntiAlpha()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theAntiAlphaInelastic);
}

void G4AntiBarionBuilder::RegisterMe(G4PhysicsBuilderInterface* aB ) {
  auto bld = dynamic_cast<G4VAntiBarionBuilder*>(aB);
  if ( bld != nullptr ) {
      theModelCollections.push_back(bld);
  } else {
      G4PhysicsBuilderInterface::RegisterMe(aB);
  }

}
