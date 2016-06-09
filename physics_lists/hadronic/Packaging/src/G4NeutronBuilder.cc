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
#include "G4NeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronBuilder::
G4NeutronBuilder(): wasActivated(false) 
{
  theNeutronElasticProcess = new G4HadronElasticProcess;
  theNeutronInelastic = new G4NeutronInelasticProcess;
  theNeutronCapture = new G4HadronCaptureProcess;
  theNeutronFission = new G4HadronFissionProcess;
  
}

G4NeutronBuilder::
~G4NeutronBuilder() 
{
  if(wasActivated)
  {
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  if(theProcMan) theProcMan->RemoveProcess(theNeutronElasticProcess);
  if(theProcMan) theProcMan->RemoveProcess(theNeutronInelastic);
  if(theProcMan) theProcMan->RemoveProcess(theNeutronCapture);
  if(theProcMan) theProcMan->RemoveProcess(theNeutronFission);
  }
  delete theNeutronElasticProcess;
  delete theNeutronInelastic;
  delete theNeutronCapture;
  delete theNeutronFission;
}

void G4NeutronBuilder::
Build()
{
  wasActivated = true;
  std::vector<G4VNeutronBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theNeutronElasticProcess);
    (*i)->Build(theNeutronInelastic);
    (*i)->Build(theNeutronCapture);
    (*i)->Build(theNeutronFission);
  }
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theNeutronElasticProcess);
  theProcMan->AddDiscreteProcess(theNeutronInelastic);
  theProcMan->AddDiscreteProcess(theNeutronCapture);
  theProcMan->AddDiscreteProcess(theNeutronFission);
}
// 2002 by J.P. Wellisch
