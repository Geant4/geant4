// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransportationManager.cc,v 1.3 1999-10-29 16:40:31 japost Exp $
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


G4TransportationManager  G4TransportationManager::fTransportationManager;

G4TransportationManager::G4TransportationManager() 
{ 
  fNavigatorForTracking= new G4Navigator() ;
  fFieldManager=         new G4FieldManager() ;
  fPropagatorInField=    new G4PropagatorInField( fNavigatorForTracking,
                                                    fFieldManager);
} 


G4TransportationManager::~G4TransportationManager()
{
  delete fNavigatorForTracking; 
  delete fPropagatorInField;
  delete fFieldManager; 
}
