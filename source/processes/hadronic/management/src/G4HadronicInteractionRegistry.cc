#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicInteraction.hh"

G4HadronicInteractionRegistry G4HadronicInteractionRegistry::
theRegistry;

void G4HadronicInteractionRegistry::
RegisterMe(G4HadronicInteraction * aModel)
{
  theRegistry.AddModel(aModel);
}

G4HadronicInteractionRegistry::~G4HadronicInteractionRegistry()
{
  for(G4int i=0; i<nModels; i++)
  {
    delete allModels[i];
  }
}

void G4HadronicInteractionRegistry::
AddModel(G4HadronicInteraction * aModel)
{
  nModels++;
  allModels.resize(nModels);
  allModels(nModels-1) = aModel;
}
