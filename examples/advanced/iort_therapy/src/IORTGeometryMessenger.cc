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
// This is the *BASIC* version of Collimator60, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTGeometryMessenger.hh"
#include "IORTGeometryController.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

IORTGeometryMessenger::IORTGeometryMessenger(IORTGeometryController* controller)
  :iortGeometryController(controller)

{
  changeTheGeometryDir = new G4UIdirectory("/geometrySetup/");
  changeTheGeometryDir -> SetGuidance("Geometry setup");

  changeTheGeometryCmd = new G4UIcmdWithAString("/geometrySetup/selectGeometry",this);
  changeTheGeometryCmd -> SetGuidance("Select the geometry you wish to use");
  changeTheGeometryCmd -> SetParameterName("Geometry",false);
//  changeTheGeometryCmd -> AvailableForStates(G4State_PreInit);
  changeTheGeometryCmd -> AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

    updateCmd = new G4UIcmdWithoutParameter("/geometrySetup/update",this);
    updateCmd->SetGuidance("Update collimator geometry.");
    updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    updateCmd->SetGuidance("if you changed geometrical value(s).");
    updateCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);
}


IORTGeometryMessenger::~IORTGeometryMessenger()
{ 
  delete changeTheGeometryDir;
  delete changeTheGeometryCmd;
  delete updateCmd;
}

void IORTGeometryMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == changeTheGeometryCmd )
    { iortGeometryController -> SetGeometry (newValue);}
 //   { iortGeometryController -> registerGeometry(G4VUserDetectorConstruction* detector);}

  else if (command == updateCmd)
  {
      iortGeometryController -> UpdateGeometry();
  }

}
