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
// Class G4TransportationManager implementation
//
// Created : J.Apostolakis, 1997
// Reviewed: G.Cosmo, 2006
// --------------------------------------------------------------------

#include "G4TransportationManager.hh"

#include <algorithm>

#include "G4GeometryMessenger.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// Initialise the static instance of the singleton
//
G4ThreadLocal G4TransportationManager*
G4TransportationManager::fTransportationManager = nullptr;

// The first registered navigator -- expect this to be the master thread's navigator
//  If it has an external sub-navigator, it will be cloned for each worker thread.
G4Navigator* G4TransportationManager::fFirstTrackingNavigator= nullptr;

// ----------------------------------------------------------------------------
// Constructor
//
G4TransportationManager::G4TransportationManager() 
{ 
  if (fTransportationManager)
  {
    G4Exception("G4TransportationManager::G4TransportationManager()",
                "GeomNav0002", FatalException,
                "Only ONE instance of G4TransportationManager is allowed!");
  }

  // Create the navigator for tracking and activate it; add to collections
  //
  G4Navigator* trackingNavigator= nullptr; 
  if( fFirstTrackingNavigator && fFirstTrackingNavigator->GetExternalNavigation() )
  {
     trackingNavigator = fFirstTrackingNavigator->Clone();
  }
  else
  {
     trackingNavigator = new G4Navigator();
     if( fFirstTrackingNavigator == nullptr ) 
        fFirstTrackingNavigator = trackingNavigator;
  }
  trackingNavigator->Activate(true);
  fNavigators.push_back(trackingNavigator);
  fActiveNavigators.push_back(trackingNavigator);
  fWorlds.push_back(trackingNavigator->GetWorldVolume()); // NULL registered

  fGeomMessenger    = new G4GeometryMessenger(this);
  fFieldManager     = new G4FieldManager(); // deleted by G4FieldManagerStore
  fPropagatorInField= new G4PropagatorInField(trackingNavigator,fFieldManager);
  fSafetyHelper     = new G4SafetyHelper();
} 

// ----------------------------------------------------------------------------
// Destructor
//
G4TransportationManager::~G4TransportationManager()
{
  delete fSafetyHelper;
  delete fPropagatorInField;
  delete fGeomMessenger;
  ClearNavigators();
  fTransportationManager = nullptr; 
}

// ----------------------------------------------------------------------------
// GetTransportationManager()
//
// Retrieve the static instance of the singleton and create it if not existing
//
G4TransportationManager* G4TransportationManager::GetTransportationManager()   
{
   if (fTransportationManager == nullptr)
   {
     fTransportationManager = new G4TransportationManager;
   }   
   return fTransportationManager;
}

// ----------------------------------------------------------------------------
// GetInstanceIfExist()
//
// Retrieve the static instance pointer of the singleton
//
G4TransportationManager* G4TransportationManager::GetInstanceIfExist()
{
   return fTransportationManager;
}

// ----------------------------------------------------------------------------
// SetFieldManager()
//
// Set the associated field manager.
//
void G4TransportationManager::SetFieldManager(G4FieldManager* newFieldManager)
{
   fFieldManager = newFieldManager; 

   // Message the PropagatorInField, 
   // which also maintains this information (to be reviewed)
   //
   if( fPropagatorInField )
   {
      fPropagatorInField -> SetDetectorFieldManager( newFieldManager );
   }
}

// ----------------------------------------------------------------------------
// SetNavigatorForTracking()
//
// Set the active navigator for tracking, always
// the first in the collection of registered navigators.
//
void G4TransportationManager::SetNavigatorForTracking(G4Navigator* newNavigator)
{
   fNavigators[0] = newNavigator;
   fActiveNavigators[0] = newNavigator;
   fPropagatorInField->SetNavigatorForPropagating(newNavigator);
}

// ----------------------------------------------------------------------------
// ClearNavigators()
//
// Clear collection of navigators and delete allocated objects.
// Called only by the class destructor.
//
void G4TransportationManager::ClearNavigators()
{
   for (auto pNav=fNavigators.cbegin(); pNav!=fNavigators.cend(); ++pNav)
   {
     delete *pNav;
   }
   fNavigators.clear();
   fActiveNavigators.clear();
   fWorlds.clear();
}

// ----------------------------------------------------------------------------
// GetParallelWorld()
//
// Provided the name of a world volume, returns the associated world pointer.
// If not existing, create (allocate) and register it in the collection.
//
G4VPhysicalVolume*
G4TransportationManager::GetParallelWorld( const G4String& worldName )
{
   G4VPhysicalVolume* wPV = IsWorldExisting(worldName);
   if (wPV == nullptr)
   {
     wPV = GetNavigatorForTracking()->GetWorldVolume();
     G4LogicalVolume* wLV = wPV->GetLogicalVolume();
     wLV = new G4LogicalVolume(wLV->GetSolid(), nullptr,
                               worldName);
     wPV = new G4PVPlacement (wPV->GetRotation(),
                              wPV->GetTranslation(),
                              wLV, worldName, nullptr, false, 0);
     RegisterWorld(wPV);
   }
   return wPV;
}

// ----------------------------------------------------------------------------
// GetNavigator()
//
// Provided the name of a world volume, returns the associated navigator.
// If not existing, create it and register it in the collection, throw an
// exception if the associated parallel world does not exist.
//
G4Navigator* G4TransportationManager::GetNavigator( const G4String& worldName )
{
   // If already existing, return the stored pointer to the navigator
   //
   for (auto pNav=fNavigators.cbegin(); pNav!=fNavigators.cend(); ++pNav)
   {
      if ((*pNav)->GetWorldVolume()->GetName() == worldName) { return *pNav; }
   }

   // Check if world of that name already exists,
   // create a navigator and register it
   //
   G4Navigator* aNavigator = nullptr;
   G4VPhysicalVolume* aWorld = IsWorldExisting(worldName);
   if(aWorld != nullptr)
   {
      aNavigator = new G4Navigator();
      aNavigator->SetWorldVolume(aWorld);
      fNavigators.push_back(aNavigator);
   }
   else
   {
      G4String message
         = "World volume with name -" + worldName
         + "- does not exist. Create it first by GetParallelWorld() method!";      
      G4Exception("G4TransportationManager::GetNavigator(name)",
                  "GeomNav0002", FatalException, message);
   }

   return aNavigator;
}

// ----------------------------------------------------------------------------
// GetNavigator()
//
// Provided a pointer to a world volume, returns the associated navigator.
// Create it in case not existing and add it to the collection.
// If world volume not existing, issue an exception.
//
G4Navigator* G4TransportationManager::GetNavigator( G4VPhysicalVolume* aWorld )
{
   for (auto pNav=fNavigators.cbegin(); pNav!=fNavigators.cend(); ++pNav)
   {
     if ((*pNav)->GetWorldVolume() == aWorld) { return *pNav; }
   }
   G4Navigator* aNavigator = nullptr;
   auto pWorld = std::find(fWorlds.cbegin(), fWorlds.cend(), aWorld);
   if (pWorld != fWorlds.cend())
   {
      aNavigator = new G4Navigator();
      aNavigator->SetWorldVolume(aWorld);
      fNavigators.push_back(aNavigator);
   }
   else
   {
      G4String message
         = "World volume with name -" + aWorld->GetName()
         + "- does not exist. Create it first by GetParallelWorld() method!";
      G4Exception("G4TransportationManager::GetNavigator(pointer)",
                  "GeomNav0002", FatalException, message);
   }

   return aNavigator;
}

// ----------------------------------------------------------------------------
// DeRegisterNavigator()
//
// Provided a pointer to an already allocated navigator object, removes the
// associated entry in the navigators collection (remove pair) but does not
// delete the actual pointed object, which is still owned by the caller.
// The navigator for tracking -cannot- be deregistered.
//
void G4TransportationManager::DeRegisterNavigator( G4Navigator* aNavigator )
{
   if (aNavigator == fNavigators[0])
   {
      G4Exception("G4TransportationManager::DeRegisterNavigator()",
                  "GeomNav0003", FatalException,
                  "The navigator for tracking CANNOT be deregistered!");
   }
   auto pNav = std::find(fNavigators.cbegin(), fNavigators.cend(), aNavigator);
   if (pNav != fNavigators.cend())
   {
      // Deregister associated world volume
      //
      DeRegisterWorld((*pNav)->GetWorldVolume());

      // Deregister the navigator
      //
      fNavigators.erase(pNav);
   }
   else
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";      
      G4Exception("G4TransportationManager::DeRegisterNavigator()",
                  "GeomNav1002", JustWarning, message);
   }
}

// ----------------------------------------------------------------------------
// ActivateNavigator()
//
// Provided a pointer to an already allocated navigator object, set to 'true'
// the associated activation flag for the navigator in the collection.
// If the provided navigator is not already registered, issue a warning
// Return the index of the activated navigator. This index should be used for
// ComputeStep() method of G4PathFinder.
//
G4int G4TransportationManager::ActivateNavigator( G4Navigator* aNavigator )
{
   auto pNav = std::find(fNavigators.cbegin(), fNavigators.cend(), aNavigator);
   if (pNav == fNavigators.cend())
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";      
      G4Exception("G4TransportationManager::ActivateNavigator()",
                  "GeomNav1002", FatalException, message);
      return -1;
   }

   aNavigator->Activate(true);
   G4int id = 0;
   for(auto pActiveNav=fActiveNavigators.cbegin();
       pActiveNav!=fActiveNavigators.cend(); ++pActiveNav)
   {
      if (*pActiveNav == aNavigator)  { return id; }
      ++id;
   }
   
   fActiveNavigators.push_back(aNavigator);
   return id;
}

// ----------------------------------------------------------------------------
// DeActivateNavigator()
//
// Provided a pointer to an already allocated navigator object, set to 'false'
// the associated activation flag in the navigators collection.
// If the provided navigator is not already registered, issue a warning.
//
void G4TransportationManager::DeActivateNavigator( G4Navigator* aNavigator )
{
   auto pNav = std::find(fNavigators.cbegin(), fNavigators.cend(), aNavigator);
   if (pNav != fNavigators.cend())
   {
      (*pNav)->Activate(false);
   }
   else
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";
      G4Exception("G4TransportationManager::DeActivateNavigator()",
                  "GeomNav1002", JustWarning, message);
   }

   auto pActiveNav = std::find(fActiveNavigators.cbegin(),
                               fActiveNavigators.cend(), aNavigator);
   if (pActiveNav != fActiveNavigators.cend())
   {
      fActiveNavigators.erase(pActiveNav);
   }
}

// ----------------------------------------------------------------------------
// InactivateAll()
//
// Inactivate all the navigators except for the tracking one, and clear the
// store of active navigators.
//
void G4TransportationManager::InactivateAll( )
{
   for (auto pNav=fActiveNavigators.cbegin();
             pNav!=fActiveNavigators.cend(); ++pNav)
   {
      (*pNav)->Activate(false);
   }
   fActiveNavigators.clear();

   // Restore status for the navigator for tracking
   //
   fNavigators[0]->Activate(true);
   fActiveNavigators.push_back(fNavigators[0]);
}

// ----------------------------------------------------------------------------
// IsWorldExisting()
//
// Verify existance or not of an istance of the world volume with
// same name in the collection. Return the world pointer if existing.
//
G4VPhysicalVolume*
G4TransportationManager::IsWorldExisting ( const G4String& name )
{
   auto pWorld = fWorlds.begin();
   if ( *pWorld==nullptr )  { *pWorld=fNavigators[0]->GetWorldVolume(); }

   for (auto cpWorld=fWorlds.cbegin(); cpWorld!=fWorlds.cend(); ++cpWorld)
   {
      if ((*cpWorld)->GetName() == name ) { return *cpWorld; }
   }
   return 0;
}

// ----------------------------------------------------------------------------
// RegisterWorld()
//
// Provided a pointer to an already allocated world object, check and add the
// associated entry in the worlds collection. Return 'true' if registration
// succeeds and the new entry is created.
//
G4bool G4TransportationManager::RegisterWorld( G4VPhysicalVolume* aWorld )
{
   G4bool done = false;

   auto pWorld = std::find(fWorlds.cbegin(), fWorlds.cend(), aWorld);
   if (pWorld == fWorlds.cend())
   {
     fWorlds.push_back(aWorld);
     done = true;
   }
   return done;
}

// ----------------------------------------------------------------------------
// DeRegisterWorld()
//
// Provided a pointer to an already allocated world object, removes the
// associated entry in the worlds collection but does not delete the actual
// pointed object, which is still owned by the caller.
//
void G4TransportationManager::DeRegisterWorld( G4VPhysicalVolume* aWorld )
{
   auto pWorld = std::find(fWorlds.cbegin(), fWorlds.cend(), aWorld);
   if (pWorld != fWorlds.cend())
   {
      fWorlds.erase(pWorld);
   }
   else
   {
     G4String message
       = "World volume -" + aWorld->GetName() + "- not found in memory!";      
     G4Exception("G4TransportationManager::DeRegisterWorld()",
                 "GeomNav1002", JustWarning, message);
   }
}

// ----------------------------------------------------------------------------
// ClearParallelWorlds()
//
// Clear collection of navigators and delete allocated objects associated with
// parallel worlds.
// Called only by the RunManager when the entire geometry is rebuilt from
// scratch.
//
void G4TransportationManager::ClearParallelWorlds()
{
   auto pNav = fNavigators.cbegin();
   G4Navigator* trackingNavigator = *pNav;
   for (pNav=fNavigators.cbegin(); pNav!=fNavigators.cend(); ++pNav)
   {
     if (*pNav != trackingNavigator)  { delete *pNav; }
   }
   fNavigators.clear();
   fActiveNavigators.clear();
   fWorlds.clear();

   fNavigators.push_back(trackingNavigator);
   fActiveNavigators.push_back(trackingNavigator);
   fWorlds.push_back(0); // NULL registered
}

// ----------------------------------------------------------------------------
// GetFirstTrackingNavigator()
//
// Get pointer to the first tracking Navigator created
// 
G4Navigator* G4TransportationManager::GetFirstTrackingNavigator()
{
  return fFirstTrackingNavigator;
}

// ----------------------------------------------------------------------------
// GetFirstTrackingNavigator()
//
// Get pointer to the first tracking Navigator created

void G4TransportationManager::SetFirstTrackingNavigator(G4Navigator *nav)
{
  fFirstTrackingNavigator= nav;
}
