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
// $Id: G4ITTransportationManager.hh 100802 2016-11-02 14:55:27Z gcosmo $
//
/// \brief {Duplicated version of G4TransportationManager.
///         This class just contains the pointer to the navigator object of the
///         simulation.}
//
// -------------------------------------------------------------------
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITTRANSPORTATIONMANAGER_HH
#define G4ITTRANSPORTATIONMANAGER_HH

#include "globals.hh"
#include "G4ITNavigator.hh"
#include <vector>

class G4ITSafetyHelper;
class G4VPhysicalVolume;

class G4ITTransportationManager
{
public:
    static void DeleteInstance();
    static G4ITTransportationManager* GetTransportationManager();

    G4ITNavigator* GetNavigatorForTracking() const;

    inline void SetWorldForTracking(G4VPhysicalVolume* theWorld);
      // Set the world volume for tracking
      // This method is to be invoked by G4RunManagerKernel.

    inline size_t GetNoActiveNavigators() const;
    inline std::vector<G4ITNavigator*>::iterator GetActiveNavigatorsIterator();
      // Return an iterator to the list of active navigators

    inline size_t GetNoWorlds() const;
    inline std::vector<G4VPhysicalVolume*>::iterator GetWorldsIterator();
      // Return an iterator to the list of registered worlds

    inline G4ITSafetyHelper* GetSafetyHelper() const;

    G4VPhysicalVolume* GetParallelWorld ( const G4String& worldName );
      // Return an exact copy of the tracking world volume. If already
      // existing just return the pointer

    G4VPhysicalVolume* IsWorldExisting ( const G4String& worldName );
      // Verify existance or not of an istance of the world volume with
      // same name in the collection

    G4ITNavigator* GetNavigator ( const G4String& worldName );
    G4ITNavigator* GetNavigator ( G4VPhysicalVolume* aWorld );
      // Return a navigator associated to either the world volume name
      // or the pointer to world physical volume. If not existing already
      // create it and register it in the collection

    G4bool RegisterWorld( G4VPhysicalVolume* aWorld );
    void DeRegisterNavigator( G4ITNavigator* aNavigator );
    G4int  ActivateNavigator( G4ITNavigator* aNavigator );
    void DeActivateNavigator( G4ITNavigator* aNavigator );
    void InactivateAll();
      // Methods for handling navigators. Navigator for tracking is always the
      // first, i.e. position 0 in the collection and cannot be de-registered

private:

   void ClearNavigators();
     // Clear collection of navigators and delete allocated objects
   void DeRegisterWorld( G4VPhysicalVolume* aWorld );
     // Register/de-register an already allocated world volume.
     // The pointed object is not deleted.

private:
    G4ITTransportationManager();
    ~G4ITTransportationManager();
    static G4ThreadLocal G4ITTransportationManager* fpInstance;
    void Initialize();
    //G4ITNavigator* fpNavigator;
    G4ITSafetyHelper* fpSafetyHelper;
   std::vector<G4ITNavigator*> fNavigators;
     // The collection of all navigators registered
   std::vector<G4ITNavigator*> fActiveNavigators;
     // The collection of only active navigators
   std::vector<G4VPhysicalVolume*> fWorlds;
     // The collection of worlds associated to the registered navigators

};

#include "G4ITTransportationManager.icc"

#endif // G4ITTRANSPORTATIONMANAGER_HH
