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
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#include "G4NeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronBuilder::
G4NeutronBuilder(G4bool fissionFlag): wasActivated(false), isFissionActivated(fissionFlag) 
{
  theNeutronInelastic = new G4NeutronInelasticProcess;
  theNeutronCapture = new G4HadronCaptureProcess;
  if ( isFissionActivated ) {
    theNeutronFission = new G4HadronFissionProcess;
  } else {
    theNeutronFission = 0;
  } 
}

void G4NeutronBuilder::
Build()
{
  wasActivated = true;
  std::vector<G4VNeutronBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theNeutronInelastic);
    (*i)->Build(theNeutronCapture);
    if ( isFissionActivated ) (*i)->Build(theNeutronFission);
  }
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theNeutronInelastic);
  theProcMan->AddDiscreteProcess(theNeutronCapture);
  if ( isFissionActivated ) theProcMan->AddDiscreteProcess(theNeutronFission);
}

void G4NeutronBuilder::RegisterMe(G4PhysicsBuilderInterface* aB) {
  auto bld = dynamic_cast<G4VNeutronBuilder*>(aB);
  if ( bld != nullptr ) {
      theModelCollections.push_back(bld);
  } else {
      G4PhysicsBuilderInterface::RegisterMe(aB);
  }
}
// 2002 by J.P. Wellisch
