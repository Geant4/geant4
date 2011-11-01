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
// $Id: AXPETDetectorMessenger.cc,v 1.2 2010-11-16 13:36:00 tnikitin Exp $
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "AXPETDetectorMessenger.hh"

#include "AXPETDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

#include "G4ios.hh"

AXPETDetectorMessenger::AXPETDetectorMessenger(AXPETDetectorConstruction * myDC)
  : myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selDetCmd = new G4UIcmdWithAString("/mydet/SelectDetector",this);
  selDetCmd->SetGuidance("Select Detector Setup.");
  selDetCmd->SetGuidance("  Choice : Detector type ");
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("Tubs");

 selDetCmd->SetCandidates("Trap Trd Tet Sphere HalfSphere HollowSphere HalfHollowSphere  Ring Shell Orb Box Cons manyCons Tubs CutTubs Hype Torus Para Paraboloid Polycone PolyconeGen PolyconeGenComplex Polyhedra PolyhedraGen PolyhedraGenComplex BREPBox Trd b1Ib2 b1Ub2 b1Sb2 b1Ub1 b1Ib1 b1Sb1 Ellipsoid EllipticalCone EllipticalTube Tet GenericTrap TwistedBox TwistedTrd TwistedTrap TwistedTubs TessellatedSolid ExtrudedSolid");
  selDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rotXCmd = new G4UIcmdWithADouble("/mydet/RotateX",this);
  rotXCmd->SetGuidance("Rotation in X direction ");
  rotXCmd->SetGuidance("  Choice : deg");
  rotXCmd->SetParameterName("choice",true);
  rotXCmd->SetDefaultValue(0.0);
  rotXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rotYCmd = new G4UIcmdWithADouble("/mydet/RotateY",this);
  rotYCmd->SetGuidance("Rotation in Y direction ");
  rotYCmd->SetGuidance("  Choice : deg");
  rotYCmd->SetParameterName("choice",true);
  rotYCmd->SetDefaultValue(0.0);
  rotYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  rotZCmd = new G4UIcmdWithADouble("/mydet/RotateZ",this);
  rotZCmd->SetGuidance("Rotation in Z direction ");
  rotZCmd->SetGuidance("  Choice : deg");
  rotZCmd->SetParameterName("choice",true);
  rotZCmd->SetDefaultValue(0.0);
  rotZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbortCmd = new G4UIcmdWithAnInteger("/mydet/AbortRun",this);
  AbortCmd->SetGuidance("Abortion of Run instead of Warnings  ");
  AbortCmd->SetGuidance("  Choice : 0/1");
  AbortCmd->SetParameterName("choice",true);
  AbortCmd->SetDefaultValue(0);
  AbortCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // myDetector->SelectDetector(defParam="Tubs");
   myDetector->SetDetectorName(defParam="Tubs");
   myDetector->Construct();




}
AXPETDetectorMessenger::~AXPETDetectorMessenger() {
  delete rotXCmd;
  delete rotYCmd;
  delete rotZCmd;
  delete AbortCmd;
  delete selDetCmd;
}
void AXPETDetectorMessenger::SetNewValue(G4UIcommand * command,
                                         G4String newValues)
{
  if( command == selDetCmd )
  {
    //myDetector->SelectDetector(newValues);
    myDetector->SetDetectorName(newValues);
    myDetector->Construct();
    myDetector->SwitchDetector();
  }
  if( command == rotXCmd)
  {
    myDetector->SetRotationInX(rotXCmd->GetNewDoubleValue(newValues));
  }
   if( command == rotYCmd)
  {
    myDetector->SetRotationInY(rotYCmd->GetNewDoubleValue(newValues));
  }
   if( command == rotZCmd)
  {
    myDetector->SetRotationInZ(rotZCmd->GetNewDoubleValue(newValues));
  }

 if( command == AbortCmd)
  {
    myDetector->SetAbortAction(AbortCmd->GetNewIntValue(newValues));
  }

  return;
}
