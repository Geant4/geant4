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
// $Id: G4TransportationManager.cc,v 1.16 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4TransportationManager 
//
// Created : J.Apostolakis, 1997
// Reviewed: G.Cosmo, 2006
//  10.04.07 V.Ivanchenko  Use unique G4SafetyHelper
//
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
G4TransportationManager* G4TransportationManager::fTransportationManager=0;

// ----------------------------------------------------------------------------
// Constructor
//
G4TransportationManager::G4TransportationManager() 
{ 
  if (fTransportationManager)
  {
    G4cerr << "Only ONE instance of G4TransportationManager is allowed!"
           << G4endl;
    G4Exception("G4TransportationManager::G4TransportationManager()",
                "InvalidSetup", FatalException,
                "Only ONE instance of G4TransportationManager is allowed!");
  }

  // Create the navigator for tracking and activate it; add to collections
  //
  G4Navigator* trackingNavigator = new G4Navigator();
  trackingNavigator->Activate(true);
  fNavigators.push_back(trackingNavigator);
  fActiveNavigators.push_back(trackingNavigator);
  fWorlds.push_back(trackingNavigator->GetWorldVolume()); // NULL registered

  fGeomMessenger    = new G4GeometryMessenger(this);
  fFieldManager     = new G4FieldManager();
  fPropagatorInField= new G4PropagatorInField(trackingNavigator,fFieldManager);
  fSafetyHelper     = new G4SafetyHelper();
} 

// ----------------------------------------------------------------------------
// Destructor
//
G4TransportationManager::~G4TransportationManager()
{
  delete fFieldManager; 
  delete fPropagatorInField;
  ClearNavigators(); 
  delete fGeomMessenger;
  delete fSafetyHelper;
}

// ----------------------------------------------------------------------------
// GetTransportationManager()
//
// Retrieve the static instance of the singleton
//
G4TransportationManager* G4TransportationManager::GetTransportationManager()
{
   static G4TransportationManager theInstance;
   if (!fTransportationManager)
     fTransportationManager = &theInstance;
   
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
// ClearNavigators()
//
// Clear collection of navigators and delete allocated objects.
// Called only by the class destructor.
//
void G4TransportationManager::ClearNavigators()
{
   std::vector<G4Navigator*>::iterator pNav;
   for (pNav=fNavigators.begin(); pNav!=fNavigators.end(); pNav++)
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
   if (!wPV)
   {
     wPV = GetNavigatorForTracking()->GetWorldVolume();
     G4LogicalVolume* wLV = wPV->GetLogicalVolume();
     wLV = new G4LogicalVolume(wLV->GetSolid(), 0,
                               worldName);
     wPV = new G4PVPlacement (wPV->GetRotation(),
                              wPV->GetTranslation(),
                              wLV, worldName, 0, false, 0);
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
   std::vector<G4Navigator*>::iterator pNav;
   for (pNav=fNavigators.begin(); pNav!=fNavigators.end(); pNav++)
   {
      if ((*pNav)->GetWorldVolume()->GetName() == worldName) { return *pNav; }
   }

   // Check if world of that name already exists,
   // create a navigator and register it
   //
   G4Navigator* aNavigator = 0;
   G4VPhysicalVolume* aWorld = IsWorldExisting(worldName);
   if(aWorld)
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
                  "InvalidSetup", FatalException, message);
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
   std::vector<G4Navigator*>::iterator pNav;
   for (pNav=fNavigators.begin(); pNav!=fNavigators.end(); pNav++)
   {
     if ((*pNav)->GetWorldVolume() == aWorld) { return *pNav; }
   }
   G4Navigator* aNavigator = 0;
   std::vector<G4VPhysicalVolume*>::iterator pWorld =
     std::find(fWorlds.begin(), fWorlds.end(), aWorld);
   if (pWorld != fWorlds.end())
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
                  "InvalidSetup", FatalException, message);
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
                  "InvalidCall", FatalException,
                  "The navigator for tracking CANNOT be deregistered!");
   }
   std::vector<G4Navigator*>::iterator pNav =
     std::find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav != fNavigators.end())
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
                  "NoEffect", JustWarning, message);
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
   std::vector<G4Navigator*>::iterator pNav =
     std::find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav == fNavigators.end())
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";      
      G4Exception("G4TransportationManager::ActivateNavigator()",
                  "NoEffect", JustWarning, message);
      return -1;
   }

   aNavigator->Activate(true);
   G4int id = 0;
   std::vector<G4Navigator*>::iterator pActiveNav;
   for(pActiveNav=fActiveNavigators.begin();
       pActiveNav!=fActiveNavigators.end(); pActiveNav++)
   {
      if (*pActiveNav == aNavigator)  { return id; }
      id++;
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
   std::vector<G4Navigator*>::iterator pNav =
     std::find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav != fNavigators.end())
   {
      (*pNav)->Activate(false);
   }
   else
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";
      G4Exception("G4TransportationManager::DeActivateNavigator()",
                  "NoEffect", JustWarning, message);
   }

   std::vector<G4Navigator*>::iterator pActiveNav =
     std::find(fActiveNavigators.begin(), fActiveNavigators.end(), aNavigator);
   if (pActiveNav != fActiveNavigators.end())
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
   std::vector<G4Navigator*>::iterator pNav;
   for (pNav=fActiveNavigators.begin(); pNav!=fActiveNavigators.end(); pNav++)
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
   std::vector<G4VPhysicalVolume*>::iterator pWorld = fWorlds.begin();
   if (*pWorld==0)  { *pWorld=fNavigators[0]->GetWorldVolume(); }

   for (pWorld=fWorlds.begin(); pWorld!=fWorlds.end(); pWorld++)
   {
      if ((*pWorld)->GetName() == name ) { return *pWorld; }
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

   std::vector<G4VPhysicalVolume*>::iterator pWorld =
     std::find(fWorlds.begin(), fWorlds.end(), aWorld);
   if (pWorld == fWorlds.end())
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
   std::vector<G4VPhysicalVolume*>::iterator pWorld =
     std::find(fWorlds.begin(), fWorlds.end(), aWorld);
   if (pWorld != fWorlds.end())
   {
      fWorlds.erase(pWorld);
   }
   else
   {
     G4String message
       = "World volume -" + aWorld->GetName() + "- not found in memory!";      
     G4Exception("G4TransportationManager::DeRegisterWorld()",
                 "InvalidSetup", FatalException, message);
   }
}
