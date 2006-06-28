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
//
// File name:     RadmonSubPhysicsListWithLabelFactory.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.cc,v 1.2 2006-06-28 13:56:49 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
