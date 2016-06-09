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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "IAEADetectorConstruction.hh"
#include "PassiveProtonBeamLine.hh"
#include "G4RunManager.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyGeometryController::HadrontherapyGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyGeometryController::~HadrontherapyGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyGeometryController::SetGeometry(G4String name)
{
    G4cout <<"Activating geometry " << name << G4endl;

    if(name == "IAEA")
    {
	registerGeometry(new IAEADetectorConstruction());
	G4cout <<"IAEA geometry activated" << G4endl;
    } 

    else if(name == "default") 
    {
	registerGeometry(new PassiveProtonBeamLine());
    } 
    else
    {
	G4cout <<"Unknown geometry: " << name << ". Geometry not changed." << G4endl;
    }
}
	
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyGeometryController::registerGeometry(G4VUserDetectorConstruction *detector)
{
	G4RunManager *runManager = G4RunManager::GetRunManager();
	runManager->SetUserInitialization(detector);
	runManager->GeometryHasBeenModified();
}

