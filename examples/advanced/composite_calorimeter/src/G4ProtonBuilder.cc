#include "G4ProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4ProtonBuilder::
G4ProtonBuilder(): wasActivated(false)  {}

G4ProtonBuilder::
~G4ProtonBuilder() 
{
  if(wasActivated)
  {
  G4ProcessManager * theProcMan = G4Proton::Proton()->GetProcessManager();
  if(theProcMan) theProcMan->RemoveProcess(&theProtonElasticProcess);
  if(theProcMan) theProcMan->RemoveProcess(&theProtonInelastic);
  }
}

void G4ProtonBuilder::
Build()
{
  wasActivated = true;
  std::vector<G4VProtonBuilder *>::iterator i;
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
