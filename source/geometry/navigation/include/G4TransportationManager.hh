//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4TransportationManager.hh 103219 2017-03-22 11:30:15Z gcosmo $
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
#include "G4SafetyHelper.hh"

#include <vector>

class G4PropagatorInField;
class G4GeometryMessenger;
class G4FieldManager;
class G4VPhysicalVolume;

class G4TransportationManager 
{
  public:  // with description

     static G4TransportationManager* GetTransportationManager();
       // Retrieve the static instance
     static G4TransportationManager* GetInstanceIfExist();
       // Retrieve singleton instance pointer.

     inline G4PropagatorInField* GetPropagatorInField() const;
     inline void SetPropagatorInField( G4PropagatorInField* newFieldPropagator );
     inline G4FieldManager* GetFieldManager() const;
     void SetFieldManager( G4FieldManager* newFieldManager );
       // Accessors for field handling

     inline G4Navigator* GetNavigatorForTracking() const;
     void SetNavigatorForTracking( G4Navigator* newNavigator );
       // Accessors for the navigator for tracking

     inline void SetWorldForTracking(G4VPhysicalVolume* theWorld);
       // Set the world volume for tracking
       // This method is to be invoked by G4RunManagerKernel.

     inline size_t GetNoActiveNavigators() const;
     inline std::vector<G4Navigator*>::iterator GetActiveNavigatorsIterator();
       // Return an iterator to the list of active navigators

     inline size_t GetNoWorlds() const;
     inline std::vector<G4VPhysicalVolume*>::iterator GetWorldsIterator();
       // Return an iterator to the list of registered worlds

     inline G4SafetyHelper* GetSafetyHelper() const;
       // Return the pointer to the navigation safety helper instance

     G4VPhysicalVolume* GetParallelWorld ( const G4String& worldName );
       // Return an exact copy of the tracking world volume. If already
       // existing just return the pointer

     G4VPhysicalVolume* IsWorldExisting ( const G4String& worldName );
       // Verify existance or not of an istance of the world volume with
       // same name in the collection

     G4Navigator* GetNavigator ( const G4String& worldName );
     G4Navigator* GetNavigator ( G4VPhysicalVolume* aWorld );
       // Return a navigator associated to either the world volume name
       // or the pointer to world physical volume. If not existing already
       // create it and register it in the collection

     G4bool RegisterWorld( G4VPhysicalVolume* aWorld );
     void DeRegisterNavigator( G4Navigator* aNavigator );
     G4int  ActivateNavigator( G4Navigator* aNavigator );
     void DeActivateNavigator( G4Navigator* aNavigator );
     void InactivateAll();
       // Methods for handling navigators. Navigator for tracking is always the
       // first, i.e. position 0 in the collection and cannot be de-registered

  public:  // without description

     void ClearParallelWorlds();
       // Clear collection of navigators and delete allocated objects
       // associated with parallel worlds. Internal method called only
       // by the RunManager when the entire geometry is rebuilt from scratch.

     ~G4TransportationManager(); 
       // Destructor

  protected:

     G4TransportationManager();
       // Singleton. Protected constructor

  private:

     void ClearNavigators();
       // Clear collection of navigators and delete allocated objects
     void DeRegisterWorld( G4VPhysicalVolume* aWorld );
       // Register/de-register an already allocated world volume.
       // The pointed object is not deleted.
 
  private:

     std::vector<G4Navigator*> fNavigators;
       // The collection of all navigators registered
     std::vector<G4Navigator*> fActiveNavigators;
       // The collection of only active navigators
     std::vector<G4VPhysicalVolume*> fWorlds;
       // The collection of worlds associated to the registered navigators

     G4PropagatorInField*    fPropagatorInField;
     G4FieldManager*         fFieldManager;
     G4GeometryMessenger*    fGeomMessenger;
     G4SafetyHelper*         fSafetyHelper;

     static G4ThreadLocal G4TransportationManager*  fTransportationManager;
};

#include "G4TransportationManager.icc"

#endif 
