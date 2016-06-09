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
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.cc,v 1.2.2.2 2006/06/29 16:13:34 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Include files
#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"
#include "RadmonDetectorLabelledEntitiesConstructorsMessenger.hh"
#include "RadmonVDetectorLabelledEntityConstructor.hh"


                                                RadmonDetectorLabelledEntitiesConstructorsFactory :: ~RadmonDetectorLabelledEntitiesConstructorsFactory()
{
 EntityConstructorsList::iterator i(entityConstructorsList.begin());
 EntityConstructorsList::iterator end(entityConstructorsList.end());
 
 RadmonDetectorLabelledEntitiesConstructorsMessenger * messenger(RadmonDetectorLabelledEntitiesConstructorsMessenger::Instance());
 
 while (i!=end)
 {
  messenger->RemoveAvailableConstructor((*i)->GetLabel());

  delete (*i);
  i++;
 }
}





RadmonVDetectorEntityConstructor *              RadmonDetectorLabelledEntitiesConstructorsFactory :: CreateEntityConstructor(const G4String & entityName)
{
 EntityConstructorsList::iterator i(entityConstructorsList.begin());
 EntityConstructorsList::iterator end(entityConstructorsList.end());
 
 while (i!=end)
 {
  if ((*i)->GetLabel()==entityName)
   return (*i)->New();

  i++;
 }
 
 return 0;
}



void                                            RadmonDetectorLabelledEntitiesConstructorsFactory :: AppendLabelledEntityConstructor(RadmonVDetectorLabelledEntityConstructor * constructor)
{
 entityConstructorsList.push_back(constructor);
 
 RadmonDetectorLabelledEntitiesConstructorsMessenger * messenger(RadmonDetectorLabelledEntitiesConstructorsMessenger::Instance());
 messenger->AddAvailableConstructor(constructor->GetLabel());
}
