#include "G4NeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronBuilder::
G4NeutronBuilder() {}

G4NeutronBuilder::
~G4NeutronBuilder() {}

void G4NeutronBuilder::
Build()
{
  G4std::vector<G4VNeutronBuilder *>::iterator i;
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
