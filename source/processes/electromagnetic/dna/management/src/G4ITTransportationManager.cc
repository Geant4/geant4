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
// $Id: G4ITTransportationManager.cc 64374 2012-10-31 16:37:23Z gcosmo $
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

G4ITTransportationManager* G4ITTransportationManager::fpInstance (0);

G4ITTransportationManager::G4ITTransportationManager()
{
    Initialize();
}

G4ITTransportationManager::~G4ITTransportationManager()
{
    if(fpNavigator) delete fpNavigator;
}

void G4ITTransportationManager::DeleteInstance()
{
    if(fpInstance)
    {
        delete fpInstance;
        fpInstance = 0;
    }
}

void G4ITTransportationManager::Initialize()
{
    fpNavigator = new G4ITNavigator();
    fpNavigator->Activate(true);
    G4Navigator* navForTracking = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4VPhysicalVolume* world = navForTracking->GetWorldVolume();
    fpNavigator->SetWorldVolume(world);
}

G4ITTransportationManager* G4ITTransportationManager::GetTransportationManager()
{
    if(fpInstance == 0) fpInstance = new G4ITTransportationManager;
    return fpInstance;
}

G4ITNavigator* G4ITTransportationManager::GetNavigatorForTracking()
{
    return fpNavigator;
}
