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
// $Id: G4TransportationManager.cc,v 1.4 2003/12/05 17:08:27 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

G4TransportationManager* G4TransportationManager::fTransportationManager=0;

G4TransportationManager::G4TransportationManager() 
{ 
  if (!fTransportationManager)
  {
    fGeomMessenger=        new G4GeometryMessenger(this);
    fNavigatorForTracking= new G4Navigator() ;
    fFieldManager=         new G4FieldManager() ;
    fPropagatorInField=    new G4PropagatorInField( fNavigatorForTracking,
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

G4TransportationManager::~G4TransportationManager()
{
  delete fGeomMessenger;
  delete fNavigatorForTracking; 
  delete fPropagatorInField;
  delete fFieldManager; 
}

G4TransportationManager* G4TransportationManager::GetTransportationManager()
{
   static G4TransportationManager theInstance;
   if (!fTransportationManager)
     fTransportationManager = &theInstance;
   
   return fTransportationManager;
}

void G4TransportationManager::SetFieldManager(G4FieldManager* newFieldManager)
{
   fFieldManager = newFieldManager; 

   // Message the PropagatorInField, 
   //     which also maintains this information (to be reviewed)
   if( fPropagatorInField )
      fPropagatorInField -> SetDetectorFieldManager( newFieldManager );
}
