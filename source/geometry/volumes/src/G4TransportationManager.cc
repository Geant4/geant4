// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransportationManager.cc,v 1.5 1999-12-15 14:50:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  G4TransportationManager 
//
//  
#include "G4TransportationManager.hh"

//  The following inclusions should be left here, as only
//    the constructor and destructor require them.
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"

// This will ensure correct order of construction and destruption of 
//  static objects.
#include "G4NavigationLevel.hh"
G4Allocator<G4NavigationLevel>     aNavigationLevelAllocator;
G4Allocator<G4NavigationLevelRep>  aNavigLevelRepAllocator;

// Ditto: correct order initialisation of class (class) data members
const G4double G4PropagatorInField::fDefault_Delta_Intersection_Val= 0.1 * mm;
const G4double G4PropagatorInField::fDefault_Delta_One_Step_Value = 0.25 * mm;


G4TransportationManager* G4TransportationManager::fTransportationManager=0;

G4TransportationManager::G4TransportationManager() 
{ 
  if (!fTransportationManager)
  {
    fNavigatorForTracking= new G4Navigator() ;
    fFieldManager=         new G4FieldManager() ;
    fPropagatorInField=    new G4PropagatorInField( fNavigatorForTracking,
                                                    fFieldManager);
  }
  else
  {
    G4Exception("Only ONE instance of G4TransportationManager is allowed!");
  }
} 


G4TransportationManager::~G4TransportationManager()
{
  delete fNavigatorForTracking; 
  delete fPropagatorInField;
  delete fFieldManager; 
}


G4TransportationManager* G4TransportationManager::GetTransportationManager()
{
   static G4TransportationManager theInstance;
   if (!fTransportationManager)
     fTransportationManager = &theInstance;
   
   return fTransportationManager;
}
