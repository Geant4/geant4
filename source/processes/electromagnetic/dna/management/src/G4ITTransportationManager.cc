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
// $Id: G4ITTransportationManager.cc 87061 2014-11-24 11:43:34Z gcosmo $
//
/// \brief {Duplicated version of G4TransportationManager.
///         This class just contains the pointer to the navigator object of the
///         simulation.}
//
// History:
// -----------
//
// -------------------------------------------------------------------

#include "G4ITTransportationManager.hh"
#include "G4TransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4ITSafetyHelper.hh"
#include "G4PVPlacement.hh"

G4ThreadLocal G4ITTransportationManager*
G4ITTransportationManager::fpInstance(0);

G4ITTransportationManager::G4ITTransportationManager()
{
  Initialize();
}

G4ITTransportationManager::~G4ITTransportationManager()
{
  ClearNavigators();
  if (fpSafetyHelper) delete fpSafetyHelper;
}

void G4ITTransportationManager::DeleteInstance()
{
  if (fpInstance)
  {
    delete fpInstance;
    fpInstance = 0;
  }
}

// ----------------------------------------------------------------------------
// ClearNavigators()
//
// Clear collection of navigators and delete allocated objects.
// Called only by the class destructor.
//
void G4ITTransportationManager::ClearNavigators()
{
  std::vector<G4ITNavigator*>::iterator pNav;
  for (pNav = fNavigators.begin(); pNav != fNavigators.end(); pNav++)
  {
    delete *pNav;
  }
  fNavigators.clear();
  fActiveNavigators.clear();
  fWorlds.clear();
}

void G4ITTransportationManager::Initialize()
{
  // Create the navigator for tracking and activate it; add to collections
  //
  G4ITNavigator* trackingNavigator = new G4ITNavigator();
  trackingNavigator->Activate(true);
  G4Navigator* navForTracking =
      G4TransportationManager::GetTransportationManager()
          ->GetNavigatorForTracking();
  G4VPhysicalVolume* world = navForTracking->GetWorldVolume();
  trackingNavigator->SetWorldVolume(world);
  fNavigators.push_back(trackingNavigator);
  fActiveNavigators.push_back(trackingNavigator);
  //fWorlds.push_back(world); // NULL registered

  size_t n_worlds = G4TransportationManager::GetTransportationManager()
      ->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator it =
      G4TransportationManager::GetTransportationManager()->GetWorldsIterator();

  for (size_t i = 0; i < n_worlds; i++, it++)
  {
    fWorlds.push_back(*it);
  }
  fpSafetyHelper = new G4ITSafetyHelper();
}

G4ITTransportationManager* G4ITTransportationManager::GetTransportationManager()
{
  if (fpInstance == 0) fpInstance = new G4ITTransportationManager;
  return fpInstance;
}

// ----------------------------------------------------------------------------
// GetParallelWorld()
//
// Provided the name of a world volume, returns the associated world pointer.
// If not existing, create (allocate) and register it in the collection.
//
G4VPhysicalVolume*
G4ITTransportationManager::GetParallelWorld(const G4String& worldName)
{
  G4VPhysicalVolume* wPV = IsWorldExisting(worldName);
  if (!wPV)
  {
    wPV = GetNavigatorForTracking()->GetWorldVolume();
    G4LogicalVolume* wLV = wPV->GetLogicalVolume();
    wLV = new G4LogicalVolume(wLV->GetSolid(), 0, worldName);
    wPV = new G4PVPlacement(wPV->GetRotation(), wPV->GetTranslation(), wLV,
                            worldName, 0, false, 0);
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
G4ITNavigator* G4ITTransportationManager::GetNavigator(const G4String& worldName)
{
  // If already existing, return the stored pointer to the navigator
  //
  std::vector<G4ITNavigator*>::iterator pNav;
  for (pNav = fNavigators.begin(); pNav != fNavigators.end(); pNav++)
  {
    if ((*pNav)->GetWorldVolume()->GetName() == worldName)
    {
      return *pNav;
    }
  }

  // Check if world of that name already exists,
  // create a navigator and register it
  //
  G4ITNavigator* aNavigator = 0;
  G4VPhysicalVolume* aWorld = IsWorldExisting(worldName);
  if (aWorld)
  {
    aNavigator = new G4ITNavigator();
    aNavigator->SetWorldVolume(aWorld);
    fNavigators.push_back(aNavigator);
  }
  else
  {
    G4String message = "World volume with name -"
        + worldName
        + "- does not exist. Create it first by GetParallelWorld() method!";
    G4Exception("G4ITTransportationManager::GetNavigator(name)", "GeomNav0002",
                FatalException, message);
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
G4ITNavigator* G4ITTransportationManager::GetNavigator(G4VPhysicalVolume* aWorld)
{
  std::vector<G4ITNavigator*>::iterator pNav;
  for (pNav = fNavigators.begin(); pNav != fNavigators.end(); pNav++)
  {
    if ((*pNav)->GetWorldVolume() == aWorld)
    {
      return *pNav;
    }
  }
  G4ITNavigator* aNavigator = 0;
  std::vector<G4VPhysicalVolume*>::iterator pWorld = std::find(fWorlds.begin(),
                                                               fWorlds.end(),
                                                               aWorld);
  if (pWorld != fWorlds.end())
  {
    aNavigator = new G4ITNavigator();
    aNavigator->SetWorldVolume(aWorld);
    fNavigators.push_back(aNavigator);
  }
  else
  {
    G4String message = "World volume with name -"
        + aWorld->GetName()
        + "- does not exist. Create it first by GetParallelWorld() method!";
    G4Exception("G4ITTransportationManager::GetNavigator(pointer)",
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
void G4ITTransportationManager::DeRegisterNavigator(G4ITNavigator* aNavigator)
{
  if (aNavigator == fNavigators[0])
  {
    G4Exception("G4ITTransportationManager::DeRegisterNavigator()",
                "GeomNav0003", FatalException,
                "The navigator for tracking CANNOT be deregistered!");
  }
  std::vector<G4ITNavigator*>::iterator pNav = std::find(fNavigators.begin(),
                                                         fNavigators.end(),
                                                         aNavigator);
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
    G4String message = "Navigator for volume -"
        + aNavigator->GetWorldVolume()->GetName() + "- not found in memory!";
    G4Exception("G4ITTransportationManager::DeRegisterNavigator()",
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
G4int G4ITTransportationManager::ActivateNavigator(G4ITNavigator* aNavigator)
{
  std::vector<G4ITNavigator*>::iterator pNav = std::find(fNavigators.begin(),
                                                         fNavigators.end(),
                                                         aNavigator);
  if (pNav == fNavigators.end())
  {
    G4String message = "Navigator for volume -"
        + aNavigator->GetWorldVolume()->GetName() + "- not found in memory!";
    G4Exception("G4ITTransportationManager::ActivateNavigator()", "GeomNav1002",
                JustWarning, message);
    return -1;
  }

  aNavigator->Activate(true);
  G4int id = 0;
  std::vector<G4ITNavigator*>::iterator pActiveNav;
  for (pActiveNav = fActiveNavigators.begin();
      pActiveNav != fActiveNavigators.end(); pActiveNav++)
  {
    if (*pActiveNav == aNavigator)
    {
      return id;
    }
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
void G4ITTransportationManager::DeActivateNavigator(G4ITNavigator* aNavigator)
{
  std::vector<G4ITNavigator*>::iterator pNav = std::find(fNavigators.begin(),
                                                         fNavigators.end(),
                                                         aNavigator);
  if (pNav != fNavigators.end())
  {
    (*pNav)->Activate(false);
  }
  else
  {
    G4String message = "Navigator for volume -"
        + aNavigator->GetWorldVolume()->GetName() + "- not found in memory!";
    G4Exception("G4ITTransportationManager::DeActivateNavigator()",
                "GeomNav1002", JustWarning, message);
  }

  std::vector<G4ITNavigator*>::iterator pActiveNav = std::find(
      fActiveNavigators.begin(), fActiveNavigators.end(), aNavigator);
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
void G4ITTransportationManager::InactivateAll()
{
  std::vector<G4ITNavigator*>::iterator pNav;
  for (pNav = fActiveNavigators.begin(); pNav != fActiveNavigators.end();
      pNav++)
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
G4ITTransportationManager::IsWorldExisting(const G4String& name)
{
  std::vector<G4VPhysicalVolume*>::iterator pWorld = fWorlds.begin();
  if (*pWorld == 0)
  {
    *pWorld = fNavigators[0]->GetWorldVolume();
  }

  for (pWorld = fWorlds.begin(); pWorld != fWorlds.end(); pWorld++)
  {
    if ((*pWorld)->GetName() == name)
    {
      return *pWorld;
    }
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
G4bool G4ITTransportationManager::RegisterWorld(G4VPhysicalVolume* aWorld)
{
  G4bool done = false;

  std::vector<G4VPhysicalVolume*>::iterator pWorld = std::find(fWorlds.begin(),
                                                               fWorlds.end(),
                                                               aWorld);
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
void G4ITTransportationManager::DeRegisterWorld(G4VPhysicalVolume* aWorld)
{
  std::vector<G4VPhysicalVolume*>::iterator pWorld = std::find(fWorlds.begin(),
                                                               fWorlds.end(),
                                                               aWorld);
  if (pWorld != fWorlds.end())
  {
    fWorlds.erase(pWorld);
  }
  else
  {
    G4String message = "World volume -" + aWorld->GetName()
                       + "- not found in memory!";
    G4Exception("G4ITTransportationManager::DeRegisterWorld()", "GeomNav1002",
                JustWarning, message);
  }
}
