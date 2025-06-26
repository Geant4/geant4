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
// G4TransportationManager
//
// Class description:
//
// A singleton class which stores the (volume) navigator used by 
// the transportation process to do the geometrical tracking.
// It also stores a pointer to the propagator used in a (magnetic) 
// field and to the field manager.
// The class instance is created before main() is called, and
// in turn creates the navigator and the rest.

// Created: John Apostolakis (CERN), 10 March 1997
// Reviewed: Gabriele Cosmo (CERN), 26 April 2006
// --------------------------------------------------------------------
#ifndef  G4TransportationManager_hh
#define  G4TransportationManager_hh 1

#include "G4Navigator.hh"
#include "G4SafetyHelper.hh"

#include <vector>

class G4PropagatorInField;
class G4GeometryMessenger;
class G4FieldManager;
class G4VPhysicalVolume;

/**
 * @brief G4TransportationManager is a singleton class which stores the
 * navigator used by the transportation process to do the geometrical tracking.
 * It also stores a pointer to the propagator used in a (magnetic) field and
 * to the field manager.
 */

class G4TransportationManager 
{
  public:

    /**
     * Retrieve the static instance.
     */
    static G4TransportationManager* GetTransportationManager();

    /**
     * Retrieve singleton instance pointer.
     */
    static G4TransportationManager* GetInstanceIfExist();

    /**
     * Accessors and modifiers for field handling.
     */
    inline G4PropagatorInField* GetPropagatorInField() const;
    inline void SetPropagatorInField(G4PropagatorInField* newFieldPropagator);
    inline G4FieldManager* GetFieldManager() const;
    void SetFieldManager( G4FieldManager* newFieldManager );

    /**
     * Accessor and modifier for the navigator for tracking.
     */
    inline G4Navigator* GetNavigatorForTracking() const;
    void SetNavigatorForTracking( G4Navigator* newNavigator );

    /**
     * Sets the world volume for tracking.
     * This method is to be invoked by G4RunManagerKernel.
     */
    inline void SetWorldForTracking(G4VPhysicalVolume* theWorld);

    /**
     * Accessors for the active navigators.
     *  @returns An iterator to the list of active navigators.
     */
    inline std::vector<G4Navigator*>::iterator GetActiveNavigatorsIterator();
    inline std::size_t GetNoActiveNavigators() const;

    /**
     * Accessors for the registered worlds.
     *  @returns An iterator to the list of registered worlds.
     */
    inline std::vector<G4VPhysicalVolume*>::iterator GetWorldsIterator();
    inline std::size_t GetNoWorlds() const;

    /**
     * Returns the pointer to the navigation safety helper instance.
     *  @returns Pointer to the navigation safety helper instance.
     */
    inline G4SafetyHelper* GetSafetyHelper() const;

    /**
     * Returns an exact copy of the tracking world volume.
     * If already existing just returns the pointer.
     *  @returns Pointer to the tracking world volume.
     */
    G4VPhysicalVolume* GetParallelWorld ( const G4String& worldName );

    /**
     * Verifies existance or not of an istance of the world volume with
     * same name in the collection.
     *  @returns Pointer to the tracking world volume.
     */
    G4VPhysicalVolume* IsWorldExisting ( const G4String& worldName );

    /**
     * Returns a navigator associated to either the world volume name
     * or associated to the pointer to the world physical volume.
     * If not existing already, creates it and registers it in the collection.
     *  @returns Pointer to a tracking navigator.
     */
    G4Navigator* GetNavigator ( const G4String& worldName );
    G4Navigator* GetNavigator ( G4VPhysicalVolume* aWorld );

    /**
     * Methods for handling navigators. Navigator for tracking is always the
     * first (i.e. position 0 in the collection) and cannot be de-registered.
     */
    G4bool RegisterWorld( G4VPhysicalVolume* aWorld );
    void DeRegisterNavigator( G4Navigator* aNavigator );
    G4int  ActivateNavigator( G4Navigator* aNavigator );
    void DeActivateNavigator( G4Navigator* aNavigator );
    void InactivateAll();

    /**
     * Accessor and modifier for the tracking navigator.
     * Retrieves/sets the first navigator pointer for the 'mass' geometry.
     * It will be used as a template for cloning the tracking navigator of
     * additional threads.
     */
    static G4Navigator* GetFirstTrackingNavigator();
    static void SetFirstTrackingNavigator(G4Navigator *nav);   

    /**
     * Clears collection of navigators and deletes the allocated objects
     * associated with parallel worlds. Internal method, called only
     * by the RunManager when the entire geometry is rebuilt from scratch.
     */
    void ClearParallelWorlds();

    /**
     * Destructor. Called internally only by G4RunManagerKernel.
     */
    ~G4TransportationManager(); 

  public:

    /** Navigator identifier. Accessed by G4CoupledTransportation. */
    static constexpr G4int kMassNavigatorId = 0;

  private:

    /**
     * Private Constructor.
     */
    G4TransportationManager();

    /**
     * Clears collection of navigators and deletes allocated objects.
     */
    void ClearNavigators();

    /**
     * De-registers an already allocated world volume.
     * The pointed object is not deleted.
     */
    void DeRegisterWorld( G4VPhysicalVolume* aWorld );
 
  private:

    /** The collection of all navigators registered. */
    std::vector<G4Navigator*> fNavigators;

    /** The collection of only active navigators. */
    std::vector<G4Navigator*> fActiveNavigators;

    /** The collection of worlds associated to the registered navigators. */
    std::vector<G4VPhysicalVolume*> fWorlds;

    G4PropagatorInField* fPropagatorInField;
    G4FieldManager*      fFieldManager;
    G4GeometryMessenger* fGeomMessenger;
    G4SafetyHelper*      fSafetyHelper;

    static G4ThreadLocal G4TransportationManager* fTransportationManager;

    static G4Navigator* fFirstTrackingNavigator;
};

#include "G4TransportationManager.icc"

#endif 
