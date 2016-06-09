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
// File name:     RadmonSubPhysicsListWithLabelFactory.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.cc,v 1.3 2006/06/29 16:20:10 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//

// Include files
#include "RadmonSubPhysicsListWithLabelFactory.hh"
#include "RadmonSubPhysicsListWithLabelMessenger.hh"
#include "RadmonVSubPhysicsListWithLabel.hh"


                                                RadmonSubPhysicsListWithLabelFactory :: ~RadmonSubPhysicsListWithLabelFactory()
{
 SubPhysicsLists::iterator i(subPhysicsLists.begin());
 SubPhysicsLists::iterator end(subPhysicsLists.end());
 
 RadmonSubPhysicsListWithLabelMessenger * messenger(RadmonSubPhysicsListWithLabelMessenger::Instance());
 
 while (i!=end)
 {
  messenger->RemoveAvailablePhysicsList((*i)->GetLabel());

  delete (*i);
  i++;
 }
}





RadmonVSubPhysicsList *                         RadmonSubPhysicsListWithLabelFactory :: CreateSubPhysicsList(const G4String & subPhysicsListName)
{
 SubPhysicsLists::iterator i(subPhysicsLists.begin());
 SubPhysicsLists::iterator end(subPhysicsLists.end());
 
 while (i!=end)
 {
  if ((*i)->GetLabel()==subPhysicsListName)
   return (*i)->New();

  i++;
 }
 
 return 0;
}



void                                            RadmonSubPhysicsListWithLabelFactory :: AppendSubPhysicsListWithLabel(RadmonVSubPhysicsListWithLabel * physicsList)
{
 subPhysicsLists.push_back(physicsList);
 
 RadmonSubPhysicsListWithLabelMessenger * messenger(RadmonSubPhysicsListWithLabelMessenger::Instance());
 messenger->AddAvailablePhysicsList(physicsList->GetLabel());
}
