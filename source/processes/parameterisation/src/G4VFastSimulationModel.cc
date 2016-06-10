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
// $Id: G4VFastSimulationModel.cc 68056 2013-03-13 14:44:48Z gcosmo $
//
//---------------------------------------------------------------
//
//  G4VFastSimulationModel.cc
//
//  Description:
//    Base class for fast simulation models.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------


#include "G4VFastSimulationModel.hh"
#include "G4FastSimulationManager.hh"

// ----------------------
// -- Simple constructor:
// ----------------------
G4VFastSimulationModel::G4VFastSimulationModel(const G4String& aName)
 : theModelName(aName) {}

// ----------------------------------------------------------------------------------------------
// -- Constructor with automatic G4FastSimulationManager constructed if needed fo given envelope:
// ----------------------------------------------------------------------------------------------
G4VFastSimulationModel::G4VFastSimulationModel(const G4String&      aName,
					       G4Envelope*     anEnvelope,
					       G4bool            IsUnique) 
 : theModelName(aName)
{
  // Retrieves the Fast Simulation Manager ou creates one if needed.
  G4FastSimulationManager* theFastSimulationManager;
  if ((theFastSimulationManager = anEnvelope->GetFastSimulationManager()) == 0) 
    theFastSimulationManager = new G4FastSimulationManager(anEnvelope,IsUnique);
  // adds this model to the Fast Simulation Manager.
  theFastSimulationManager->AddFastSimulationModel(this);
}
