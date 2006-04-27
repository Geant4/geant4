//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TransportationManager.hh,v 1.4 2006-04-27 16:30:05 gcosmo Exp $
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

// Created:  10 March 1997, J. Apostolakis
// Reviewed: 26 April 2006, G. Cosmo
// --------------------------------------------------------------------

#ifndef  G4TransportationManager_hh
#define  G4TransportationManager_hh

#include "G4Navigator.hh"

#include <vector>

class G4PropagatorInField;
class G4GeometryMessenger;
class G4FieldManager;

class G4TransportationManager 
{
  public:  // with description

     static G4TransportationManager* GetTransportationManager();
       // Retrieve the static instance

     inline G4PropagatorInField* GetPropagatorInField() const;
     inline void SetPropagatorInField( G4PropagatorInField* newFieldPropagator );
     inline G4FieldManager* GetFieldManager() const;
     void SetFieldManager( G4FieldManager* newFieldManager );
       // Accessors for field handling

     inline G4Navigator* GetNavigatorForTracking() const;
     inline void SetNavigatorForTracking( G4Navigator* newNavigator );
       // Accessors for the navigator for tracking

     inline std::vector<G4Navigator*>::iterator GetActiveNavigatorsIterator();
       // Return an iterator to the list of active navigators

     G4bool RegisterNavigator( G4Navigator* aNavigator );
     void DeRegisterNavigator( G4Navigator* aNavigator );
     G4int  ActivateNavigator( G4Navigator* aNavigator );
     void DeActivateNavigator( G4Navigator* aNavigator );
       // Methods for handling navigators. Navigator for tracking is always the
       // first, i.e. position 0 in the collection and cannot be deregistered

  protected:

     G4TransportationManager();
     ~G4TransportationManager(); 
       // Singleton. Protected constructor and destructor

  private:

     void ClearNavigators();
       // Clear collection of navigators and delete allocated objects

  private:

     std::vector<G4Navigator*> fNavigators;
       // The collection of all navigators registered
     std::vector<G4Navigator*> fActiveNavigators;
       // The collection of only active navigators

     G4PropagatorInField*    fPropagatorInField;
     G4FieldManager*         fFieldManager;
     G4GeometryMessenger*    fGeomMessenger;

     static G4TransportationManager*  fTransportationManager;
};

#include "G4TransportationManager.icc"

#endif 
