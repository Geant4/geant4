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
/// \file field/BlineTracer/src/G4BlineTracerMessenger.cc
/// \brief Implementation of the G4BlineTracerMessenger class
//
//
//
// --------------------------------------------------------------------
//
// G4BlineTracerMessenger implementation
//
// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------

#include "G4BlineTracerMessenger.hh"
#include "G4BlineTracer.hh"
#include "G4BlineEventAction.hh"

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"

///////////////////////////////////////////////////////////////////////////

G4BlineTracerMessenger::G4BlineTracerMessenger( G4BlineTracer* aBlineTool )
{
  fTheBlineTool = aBlineTool;
  fBlineToolDir = new G4UIdirectory("/vis/blineTracer/");
  fBlineToolDir->SetGuidance("Commands to trace and visualise magnetic field lines.");
  fBlineToolDir->SetGuidance("These commands work only if a magnetic-field is set");
  fBlineToolDir->SetGuidance("in the application.");

  // commands

  fBlineCmd = new G4UIcmdWithAnInteger("/vis/blineTracer/computeBline",this);
  fBlineCmd->SetGuidance("Compute magnetic field lines for visualisation.");
  fBlineCmd->SetParameterName("nb_of_lines",false);
  fBlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSetMaxTrackingStepCmd =
    new G4UIcmdWithADoubleAndUnit("/vis/blineTracer/setMaxStepLength",this); 
  fSetMaxTrackingStepCmd->SetGuidance("Set the maximum length of tracking step");
  fSetMaxTrackingStepCmd->SetGuidance("when integrating magnetic field line.");
  fSetMaxTrackingStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSetDrawColourCmd = new G4UIcmdWith3Vector("/vis/blineTracer/setColour",this);
  fSetDrawColourCmd->SetGuidance("Set the colour drawing trajectories");
  fSetDrawColourCmd->SetGuidance("and magnetic field lines.");
  fSetDrawColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSetDrawBlineCmd = new G4UIcmdWithABool("/vis/blineTracer/stockLines",this);
  fSetDrawBlineCmd->SetGuidance("If true field lines are stocked in lines");
  fSetDrawBlineCmd->SetGuidance("to be drawn.");
  fSetDrawBlineCmd->SetParameterName("StockLines",false);
  fSetDrawBlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSetDrawPointsCmd = new G4UIcmdWithABool("/vis/blineTracer/stockPoints",this);
  fSetDrawPointsCmd->SetGuidance("If true step field line points are stocked");
  fSetDrawPointsCmd->SetGuidance("in vector of points to be drawn.");
  fSetDrawPointsCmd->SetParameterName("StockPoints",false);
  fSetDrawPointsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSetPointSizeCmd = new G4UIcmdWithADouble("/vis/blineTracer/setPointSize",this);
  fSetPointSizeCmd->SetGuidance("Set the size of points for drawing.");
  fSetPointSizeCmd->SetParameterName("StepSize",false);
  fSetPointSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDrawCmd = new G4UIcmdWithoutParameter("/vis/blineTracer/show",this);
  fDrawCmd->SetGuidance("Show the stored magnetic field lines.");
  fDrawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fResetCmd =
    new G4UIcmdWithoutParameter("/vis/blineTracer/resetMaterialToBeDrawn",this);
  fResetCmd->SetGuidance("Clear the vectors of lines and points to be drawn.");
  fResetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

///////////////////////////////////////////////////////////////////////////

G4BlineTracerMessenger::~G4BlineTracerMessenger()
{
  delete fResetCmd;
  delete fDrawCmd;
  delete fSetPointSizeCmd;
  delete fSetDrawPointsCmd;
  delete fSetDrawBlineCmd;
  delete fSetDrawColourCmd;
  delete fSetMaxTrackingStepCmd;
  delete fBlineCmd;
  delete fBlineToolDir;
}  

///////////////////////////////////////////////////////////////////////////

void G4BlineTracerMessenger::SetNewValue( G4UIcommand * command,
                                          G4String newValues )
{ 
  if (command == fBlineCmd)
    fTheBlineTool->ComputeBlines(1);
  else if( command == fSetMaxTrackingStepCmd ) 
    fTheBlineTool->SetMaxTrackingStep(fSetMaxTrackingStepCmd
                                   ->GetNewDoubleValue(newValues));
  else if( command == fSetDrawBlineCmd ) 
    fTheBlineTool->GetEventAction()->SetDrawBline(fSetDrawBlineCmd
                                   ->GetNewBoolValue(newValues));
  else if( command == fSetDrawColourCmd ) 
  {
    G4ThreeVector vec=fSetDrawColourCmd->GetNew3VectorValue(newValues);
    fTheBlineTool->GetEventAction()->
                  SetDrawColour(G4Colour(vec.x(),vec.y(),vec.z()));
  }
  else if( command == fSetDrawPointsCmd ) 
    fTheBlineTool->GetEventAction()->SetDrawPoints(fSetDrawPointsCmd
                                   ->GetNewBoolValue(newValues));
  else if( command == fSetPointSizeCmd ) 
    fTheBlineTool->GetEventAction()->SetPointSize(fSetPointSizeCmd
                                   ->GetNewDoubleValue(newValues));      
  else if( command == fDrawCmd ) 
    fTheBlineTool->GetEventAction()->DrawFieldLines(.5,45.,45.);
  else if( command == fResetCmd )
    fTheBlineTool->GetEventAction()->ResetVectorObjectToBeDrawn();
}
