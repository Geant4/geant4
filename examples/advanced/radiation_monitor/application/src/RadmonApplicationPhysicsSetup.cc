//
// File name:     RadmonApplicationPhysicsSetup.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.cc,v 1.1 2005-11-07 17:53:48 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationPhysicsSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonSubPhysicsListWithLabelFactory.hh"


#define DECLARE_GENERATOR_CONSTRUCTOR(name)     subPhysicsList=new name();                                                               \
                                                if (subPhysicsList==0)                                                                   \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendSubPhysicsListWithLabel(subPhysicsList)

G4bool                                          RadmonApplicationPhysicsSetup :: CreateSubPhysicsList(RadmonSubPhysicsListWithLabelFactory * /* factory */)
{
// RadmonVSubPhysicsListWithLabel * subPhysicsList;

 return true;
}
