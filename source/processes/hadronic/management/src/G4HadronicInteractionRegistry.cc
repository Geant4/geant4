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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicInteraction.hh"

G4HadronicInteractionRegistry G4HadronicInteractionRegistry::
theRegistry;

void G4HadronicInteractionRegistry::
RegisterMe(G4HadronicInteraction * aModel)
{
  theRegistry.AddModel(aModel);
}

void G4HadronicInteractionRegistry::
RemoveMe(G4HadronicInteraction * aModel)
{
  theRegistry.allModels.erase(G4std::find(theRegistry.allModels.begin(), theRegistry.allModels.end(), aModel));
  theRegistry.nModels = theRegistry.allModels.size();
}

G4HadronicInteractionRegistry::~G4HadronicInteractionRegistry()
{
  while(allModels.size()!=0)
  {
    delete allModels.front();
  }
}

void G4HadronicInteractionRegistry::
AddModel(G4HadronicInteraction * aModel)
{
  G4bool alreadyThere = false;
  for(G4int i=0; i<nModels; i++)
  {
    if(allModels[i]==aModel)
    {
      alreadyThere = true;
      break;
    }
  }
  if(!alreadyThere)
  {
    nModels++;
    allModels.push_back(aModel);
  }
}
