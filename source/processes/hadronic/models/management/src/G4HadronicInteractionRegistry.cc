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
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicInteraction.hh"

G4HadronicInteractionRegistry & G4HadronicInteractionRegistry::theRegistry()
{
  static G4HadronicInteractionRegistry theRegistryInstance;
  return theRegistryInstance;
}

void G4HadronicInteractionRegistry::
RegisterMe(G4HadronicInteraction * aModel)
{
  theRegistry().AddModel(aModel);
}

void G4HadronicInteractionRegistry::
RemoveMe(G4HadronicInteraction * aModel)
{
  theRegistry().allModels.erase(std::find(theRegistry().allModels.begin(),
                                          theRegistry().allModels.end(), aModel));
  theRegistry().nModels = theRegistry().allModels.size();
}

G4HadronicInteractionRegistry::~G4HadronicInteractionRegistry()
{
  /*
  while(allModels.size()!=0)
  {
    delete allModels.front();
  }
  */
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
