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
// $Id: G4VFastSimulationModel.cc,v 1.5 2001-07-11 10:08:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

//
// Simple constructor.
//
G4VFastSimulationModel::
G4VFastSimulationModel(const G4String& aName)  : theModelName(aName) {}

//
// Constructor for beginners. We do all the job, no matters...
//
G4VFastSimulationModel::
G4VFastSimulationModel(const G4String& aName,
		       G4Envelope* anEnvelope,
		       G4bool IsUnique) : theModelName(aName)
{
  // Retrieves the Fast Simulation Manager ou creates one 
  // if needed.
  G4FastSimulationManager* theFastSimulationManager;
  if ((theFastSimulationManager=anEnvelope->GetFastSimulationManager()) == 0) 
    theFastSimulationManager= 
      new G4FastSimulationManager(anEnvelope,IsUnique);
  // adds this model to the Fast Simulation Manager.
  theFastSimulationManager->AddFastSimulationModel(this);
}
