// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN01PhysicsList.cc,v 1.1 1999-01-07 16:05:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN01PhysicsList.hh"
#include "G4ParticleTypes.hh"


ExN01PhysicsList::ExN01PhysicsList()
{;}

ExN01PhysicsList::~ExN01PhysicsList()
{;}

void ExN01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4Geantino::GeantinoDefinition();
}

void ExN01PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
}

void ExN01PhysicsList::SetCuts(G4double cut)
{
  // uppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                               
  SetVerboseLevel(0);                                                           
                                                                                
  // set cut values for gamma at first and for e- second and next for e+,       
  // because some processes for e+/e- need cut values for gamma                 
  SetCutValue(cut, "gamma");                                                    
  SetCutValue(cut, "e-");                                                       
  SetCutValue(cut, "e+");                                                       
                                                                                
  // set cut values for proton and anti_proton before all other hadrons         
  // because some processes for hadrons need cut values for proton/anti_proton  
  SetCutValue(cut, "proton");                                                   
  SetCutValue(cut, "anti_proton");                                              
                                                                                
  // set a cut value to all particles

  SetCutValueForOthers(cut);
                          
  SetVerboseLevel(temp);  
}

