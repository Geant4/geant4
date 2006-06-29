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
// $Id: G4BlineTracerMessenger.cc,v 1.2 2006-06-29 17:15:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  theBlineTool = aBlineTool;
  BlineToolDir = new G4UIdirectory("/vis/blineTracer/");
  BlineToolDir->SetGuidance("Commands to trace and visualise magnetic field lines.");
  BlineToolDir->SetGuidance("These commands work only if a magnetic-field is set");
  BlineToolDir->SetGuidance("in the application.");

  // commands

  BlineCmd = new G4UIcmdWithAnInteger("/vis/blineTracer/computeBline",this);
  BlineCmd->SetGuidance("Compute magnetic field lines for visualisation.");
  BlineCmd->SetParameterName("nb_of_lines",false);
  BlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetMaxTrackingStepCmd =
    new G4UIcmdWithADoubleAndUnit("/vis/blineTracer/setMaxStepLength",this); 
  SetMaxTrackingStepCmd->SetGuidance("Set the maximum length of tracking step");
  SetMaxTrackingStepCmd->SetGuidance("when integrating magnetic field line.");
  SetMaxTrackingStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawColourCmd = new G4UIcmdWith3Vector("/vis/blineTracer/setColour",this);
  SetDrawColourCmd->SetGuidance("Set the colour drawing trajectories");
  SetDrawColourCmd->SetGuidance("and magnetic field lines.");
  SetDrawColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetDrawBlineCmd = new G4UIcmdWithABool("/vis/blineTracer/stockLines",this);
  SetDrawBlineCmd->SetGuidance("If true field lines are stocked in lines");
  SetDrawBlineCmd->SetGuidance("to be drawn.");
  SetDrawBlineCmd->SetParameterName("StockLines",false);
  SetDrawBlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawPointsCmd = new G4UIcmdWithABool("/vis/blineTracer/stockPoints",this);
  SetDrawPointsCmd->SetGuidance("If true step field line points are stocked");
  SetDrawPointsCmd->SetGuidance("in vector of points to be drawn.");
  SetDrawPointsCmd->SetParameterName("StockPoints",false);
  SetDrawPointsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPointSizeCmd = new G4UIcmdWithADouble("/vis/blineTracer/setPointSize",this);
  SetPointSizeCmd->SetGuidance("Set the size of points for drawing.");
  SetPointSizeCmd->SetParameterName("StepSize",false);
  SetPointSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DrawCmd = new G4UIcmdWithoutParameter("/vis/blineTracer/show",this);
  DrawCmd->SetGuidance("Show the stored magnetic field lines.");
  DrawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ResetCmd =
    new G4UIcmdWithoutParameter("/vis/blineTracer/resetMaterialToBeDrawn",this);
  ResetCmd->SetGuidance("Clear the vectors of lines and points to be drawn.");
  ResetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

///////////////////////////////////////////////////////////////////////////

G4BlineTracerMessenger::~G4BlineTracerMessenger()
{
  delete ResetCmd;
  delete DrawCmd;
  delete SetPointSizeCmd;
  delete SetDrawPointsCmd;
  delete SetDrawBlineCmd;
  delete SetDrawColourCmd;
  delete SetMaxTrackingStepCmd;
  delete BlineCmd;
  delete BlineToolDir;
}  

///////////////////////////////////////////////////////////////////////////

void G4BlineTracerMessenger::SetNewValue( G4UIcommand * command,
                                          G4String newValues )
{ 
  if (command == BlineCmd)
    theBlineTool->ComputeBlines(1);
  else if( command == SetMaxTrackingStepCmd ) 
    theBlineTool->SetMaxTrackingStep(SetMaxTrackingStepCmd
                                   ->GetNewDoubleValue(newValues));
  else if( command == SetDrawBlineCmd ) 
    theBlineTool->GetEventAction()->SetDrawBline(SetDrawBlineCmd
                                   ->GetNewBoolValue(newValues));
  else if( command == SetDrawColourCmd ) 
  {
    G4ThreeVector vec=SetDrawColourCmd->GetNew3VectorValue(newValues);
    theBlineTool->GetEventAction()->
                  SetDrawColour(G4Colour(vec.x(),vec.y(),vec.z()));
  }
  else if( command == SetDrawPointsCmd ) 
    theBlineTool->GetEventAction()->SetDrawPoints(SetDrawPointsCmd
                                   ->GetNewBoolValue(newValues));
  else if( command == SetPointSizeCmd ) 
    theBlineTool->GetEventAction()->SetPointSize(SetPointSizeCmd
                                   ->GetNewDoubleValue(newValues));      
  else if( command == DrawCmd ) 
    theBlineTool->GetEventAction()->DrawFieldLines(.5,45.,45.);
  else if( command == ResetCmd )
    theBlineTool->GetEventAction()->ResetVectorObjectToBeDrawn();
}
