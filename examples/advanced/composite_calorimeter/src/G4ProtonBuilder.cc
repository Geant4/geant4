#include "G4ProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4ProtonBuilder::
G4ProtonBuilder() {}

G4ProtonBuilder::
~G4ProtonBuilder() {}

void G4ProtonBuilder::
Build()
{
  G4std::vector<G4VProtonBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theProtonElasticProcess);
    (*i)->Build(theProtonInelastic);
  }
  G4ProcessManager * theProcMan = G4Proton::Proton()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theProtonElasticProcess);
  theProcMan->AddDiscreteProcess(&theProtonInelastic);
}
// 2002 by J.P. Wellisch
