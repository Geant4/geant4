#include "G4NeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronBuilder::
G4NeutronBuilder(): wasActivated(false) {}

G4NeutronBuilder::
~G4NeutronBuilder() 
{
  if(wasActivated)
  {
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  if(theProcMan) theProcMan->RemoveProcess(&theNeutronElasticProcess);
  if(theProcMan) theProcMan->RemoveProcess(&theNeutronInelastic);
  if(theProcMan) theProcMan->RemoveProcess(&theNeutronCapture);
  if(theProcMan) theProcMan->RemoveProcess(&theNeutronFission);
  }
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
  theProcMan->AddDiscreteProcess(&theNeutronElasticProcess);
  theProcMan->AddDiscreteProcess(&theNeutronInelastic);
  theProcMan->AddDiscreteProcess(&theNeutronCapture);
  theProcMan->AddDiscreteProcess(&theNeutronFission);
}
// 2002 by J.P. Wellisch
