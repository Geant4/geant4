// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransportationManager.hh,v 1.5 2000-11-01 16:51:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4TransportationManager
//
// Class description:
//
// A singleton class which stores the (volume) navigator used by 
// the transportation process to do the geometrical tracking.
// It also stores a pointer to the propagator used in a (magnetic) 
// field and to the field manager.
// The class instance is created before main() is called, and
// in turn creates the navigator and the rest.

// =======================================================================
// Created:  10 March 1997, J. Apostolakis
// =======================================================================

#ifndef  G4TransportationManager_hh
#define  G4TransportationManager_hh

#include "G4Navigator.hh"

class G4PropagatorInField;
class G4FieldManager;

class G4TransportationManager 
{
  public:  // with description

     static G4TransportationManager* GetTransportationManager();

     inline G4Navigator*          GetNavigatorForTracking() const;
     inline G4PropagatorInField*  GetPropagatorInField() const;
     inline G4FieldManager*       GetFieldManager() const;

     inline void  SetNavigatorForTracking( G4Navigator*        newNavigator);
     inline void  SetPropagatorInField(  G4PropagatorInField*  newFieldPropagator);
     inline void  SetFieldManager(       G4FieldManager*       newFieldManager);

     ~G4TransportationManager(); 

  protected:

     G4TransportationManager(); 
    
  private:

     G4Navigator*            fNavigatorForTracking ;
     G4PropagatorInField*    fPropagatorInField;
     G4FieldManager*         fFieldManager;

     static G4TransportationManager*  fTransportationManager;
};

#include "G4TransportationManager.icc"

#endif 
