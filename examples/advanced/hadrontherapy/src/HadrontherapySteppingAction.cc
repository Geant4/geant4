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
// $Id: HadrontherapySteppingAction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyDetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyRunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include <iomanip.h>
#include "G4UImanager.hh"

HadrontherapySteppingAction::HadrontherapySteppingAction()
{
}

HadrontherapySteppingAction::~HadrontherapySteppingAction()
{
}

void HadrontherapySteppingAction::UserSteppingAction(const G4Step*)
{ 
}





