//
// File name:     RadmonPhysicsList.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsList.cc,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonPhysicsList.hh"

#include "RadmonVPhysicsLayout.hh"
#include "RadmonVSubPhysicsList.hh"
#include "RadmonVSubPhysicsListFactory.hh"

                                                RadmonPhysicsList :: RadmonPhysicsList(RadmonVPhysicsLayout * layout, RadmonVSubPhysicsListFactory * factory)
:
 physicsLayout(layout),
 constructorsFactory(factory)
{
}



                                                RadmonPhysicsList :: ~RadmonPhysicsList()
{
 Destruct();
 delete constructorsFactory;
}





void                                            RadmonPhysicsList :: OnLayoutChange(void)
{
 Destruct();
 
 // TO DO
 G4cout << "RadmonPhysicsList::OnLayoutChange: NOT IMPLEMENTED YET." << G4endl;
}





void                                            RadmonPhysicsList :: ConstructParticle(void)
{
 // TO DO
 G4cout << "RadmonPhysicsList::ConstructParticle: NOT IMPLEMENTED YET." << G4endl;
}



void                                            RadmonPhysicsList :: ConstructProcess(void)
{
 // TO DO
 G4cout << "RadmonPhysicsList::ConstructProcess: NOT IMPLEMENTED YET." << G4endl;
}



void                                            RadmonPhysicsList :: SetCuts(void)
{
 // TO DO
 G4cout << "RadmonPhysicsList::SetCuts: NOT IMPLEMENTED YET." << G4endl;
}





void                                            RadmonPhysicsList :: Destruct(void)
{
 SubPhysiscsLists::iterator i(subPhysiscsLists.begin());
 const SubPhysiscsLists::iterator end(subPhysiscsLists.end());
 
 while (i!=end)
 {
  delete (*i);
  
  i++;
 }
 
 subPhysiscsLists.clear();
}
