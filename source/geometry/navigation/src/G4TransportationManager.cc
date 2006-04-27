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
// $Id: G4TransportationManager.cc,v 1.6 2006-04-27 16:30:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4TransportationManager 
//
// Author: J.Apostolakis (John.Apostolakis@cern.ch), 1997
//
// --------------------------------------------------------------------

#include "G4TransportationManager.hh"

#include "G4GeometryMessenger.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"

// Initialise the static instance of the singleton
//
G4TransportationManager* G4TransportationManager::fTransportationManager=0;

// ----------------------------------------------------------------------------
// Constructor
//
G4TransportationManager::G4TransportationManager() 
{ 
  if (!fTransportationManager)
  {
    // Create the navigator for tracking and activate it; add to collections
    //
    G4Navigator* trackingNavigator = new G4Navigator();
    trackingNavigator->Activate(true);
    fNavigators.push_back(trackingNavigator);
    fActiveNavigators.push_back(trackingNavigator);

    fGeomMessenger     = new G4GeometryMessenger(this);
    fFieldManager      = new G4FieldManager();
    fPropagatorInField = new G4PropagatorInField(trackingNavigator,
                                                 fFieldManager);
  }
  else
  {
    G4cerr << "Only ONE instance of G4TransportationManager is allowed!"
           << G4endl;
    G4Exception("G4TransportationManager::G4TransportationManager()",
                "InvalidSetup", FatalException,
                "Only ONE instance of G4TransportationManager is allowed!");
  }
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
}

// ----------------------------------------------------------------------------
// GetTransportationManager
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
}

// ----------------------------------------------------------------------------
// RegisterNavigator()
//
// Provided a pointer to a newly allocated navigator object, registers an
// associated entry in the navigators collection with activation flag
// set to false. If the pointer is already registered, issues a warning
// and return 'false' as state.
//
G4bool G4TransportationManager::RegisterNavigator( G4Navigator* aNavigator )
{
   G4bool done = false;
   std::vector<G4Navigator*>::iterator pNav =
        find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav == fNavigators.end())
   {
      fNavigators.push_back(aNavigator);
      done = true;
   }
   return done;
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
        find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav != fNavigators.end())
   {
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
// Return the number of navigators currently activated including the new one.
//
G4int G4TransportationManager::ActivateNavigator( G4Navigator* aNavigator )
{
   std::vector<G4Navigator*>::iterator pNav =
        find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav == fNavigators.end())
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";      
      G4Exception("G4TransportationManager::ActivateNavigator()",
                  "NoEffect", JustWarning, message);
   }
   else
   {
      (*pNav)->Activate(true);
   }
   std::vector<G4Navigator*>::iterator pActiveNav =
        find(fActiveNavigators.begin(), fActiveNavigators.end(), aNavigator);
   if (pActiveNav == fActiveNavigators.end())
   {
      fActiveNavigators.push_back(aNavigator);
   }

   return fActiveNavigators.size();
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
        find(fNavigators.begin(), fNavigators.end(), aNavigator);
   if (pNav == fNavigators.end())
   {
      G4String message
         = "Navigator for volume -" + aNavigator->GetWorldVolume()->GetName()
         + "- not found in memory!";
      G4Exception("G4TransportationManager::DeActivateNavigator()",
                  "NoEffect", JustWarning, message);
   }
   else
   {
      (*pNav)->Activate(false);
   }
   std::vector<G4Navigator*>::iterator pActiveNav =
        find(fActiveNavigators.begin(), fActiveNavigators.end(), aNavigator);
   if (pActiveNav != fActiveNavigators.end())
   {
      fActiveNavigators.erase(pActiveNav);
   }
}
