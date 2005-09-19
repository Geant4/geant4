//
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.cc,v 1.2 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
