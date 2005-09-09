//
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"
#include "RadmonVDetectorLabelledEntityConstructor.hh"


                                                RadmonDetectorLabelledEntitiesConstructorsFactory :: ~RadmonDetectorLabelledEntitiesConstructorsFactory()
{
 EntityConstructorsList::iterator i(entityConstructorsList.begin());
 EntityConstructorsList::iterator end(entityConstructorsList.end());
 
 while (i!=end)
 {
  delete (*i);
  i++;
 }
}





RadmonVDetectorEntityConstructor *              RadmonDetectorLabelledEntitiesConstructorsFactory :: GetEntityConstructor(const G4String & entityName)
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
}
