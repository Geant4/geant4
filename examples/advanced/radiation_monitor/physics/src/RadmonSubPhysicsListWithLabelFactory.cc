//
// File name:     RadmonSubPhysicsListWithLabelFactory.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.cc,v 1.1 2005-11-07 17:52:36 capra Exp $
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
