//
// File name:     RadmonPhysicsDummyPhysicsList.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorDummyPhysicsList.cc,v 1.1 2005-09-14 12:31:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonPhysicsDummyPhysicsList.hh"

#include "G4Geantino.hh"


void                                            RadmonPhysicsDummyPhysicsList :: ConstructParticle(void)
{
 G4Geantino::GeantinoDefinition();
}



void                                            RadmonPhysicsDummyPhysicsList :: ConstructProcess(void)
{
 AddTransportation();
}



void                                            RadmonPhysicsDummyPhysicsList :: SetCuts(void)
{
 G4int verbosity(GetVerboseLevel());

 SetVerboseLevel(0);                                                           
 SetCutsWithDefault();   
 SetVerboseLevel(verbosity);  
}
