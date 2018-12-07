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
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorMessenger ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelDetectorMessenger.hh"

#include "GammaRayTelDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorMessenger::GammaRayTelDetectorMessenger(GammaRayTelDetectorConstruction * GammaRayTelDet)
  :GammaRayTelDetector(GammaRayTelDet)

{ 
  GammaRayTeldetDir = new G4UIdirectory("/payload/");
  GammaRayTeldetDir->SetGuidance("GammaRayTel payload control.");
  
  // converter material command
  
  ConverterMaterCmd = new G4UIcmdWithAString("/payload/setConvMat",this);
  ConverterMaterCmd->SetGuidance("Select Material of the Converter.");
  ConverterMaterCmd->SetParameterName("choice",false);
  ConverterMaterCmd->AvailableForStates(G4State_Idle);

  // converter thickness command
  
  ConverterThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setConvThick",this);
  ConverterThickCmd->SetGuidance("Set Thickness of the Converter");
  ConverterThickCmd->SetParameterName("Size",false);
  ConverterThickCmd->SetRange("Size>=0.");
  ConverterThickCmd->SetUnitCategory("Length");
  ConverterThickCmd->AvailableForStates(G4State_Idle);

  // tracker silicon thickness command

  SiliconThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiThick",this);
  SiliconThickCmd->SetGuidance("Set Thickness of the Silicon");
  SiliconThickCmd->SetParameterName("Size",false);
  SiliconThickCmd->SetRange("Size>=0.");
  SiliconThickCmd->SetUnitCategory("Length");
  SiliconThickCmd->AvailableForStates(G4State_Idle);

  // tracker silicon pitch command

  SiliconPitchCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiPitch",this);
  SiliconPitchCmd->SetGuidance("Set Pitch of the Silicon Strips");
  SiliconPitchCmd->SetParameterName("Size",false);
  SiliconPitchCmd->SetRange("Size>=0."); 
  SiliconPitchCmd->SetUnitCategory("Length");
  SiliconPitchCmd->AvailableForStates(G4State_Idle);
  
  // tracker silicon tile size command
  
  SiliconTileXYCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiTileXY",this);
  SiliconTileXYCmd->SetGuidance("Set XY dimensions of Si Tile");
  SiliconTileXYCmd->SetParameterName("Size",false);
  SiliconTileXYCmd->SetRange("Size>=0.");
  SiliconTileXYCmd->SetUnitCategory("Length");  
  SiliconTileXYCmd->AvailableForStates(G4State_Idle);
  
  // tracker number of silicon tiles

  NbSiTilesCmd = new G4UIcmdWithAnInteger("/payload/setNbOfSiTiles",this);
  NbSiTilesCmd->SetGuidance("Set number of Si Tiles.");
  NbSiTilesCmd->SetParameterName("NbSiTiles",false);
  NbSiTilesCmd->SetRange("NbSiTiles>0 && NbSiTiles<100");
  NbSiTilesCmd->AvailableForStates(G4State_Idle);

  // tracker number of silicon layers

  NbTKRLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfTKRLayers",this);
  NbTKRLayersCmd->SetGuidance("Set number of TKR Layers.");
  NbTKRLayersCmd->SetParameterName("NbTKRLayers",false);
  NbTKRLayersCmd->SetRange("NbTKRLayers>0 && NbTKRLayers<30");
  NbTKRLayersCmd->AvailableForStates(G4State_Idle);

  // tracker layer distance

  LayerDistanceCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setLayerDistance",this);
  LayerDistanceCmd->SetGuidance("Set distance between two layers");
  LayerDistanceCmd->SetParameterName("Size",false);
  LayerDistanceCmd->SetRange("Size>=0.");
  LayerDistanceCmd->SetUnitCategory("Length");  
  LayerDistanceCmd->AvailableForStates(G4State_Idle);

  // tracker views distance
  
  ViewsDistanceCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setViewsDistance",this);
  ViewsDistanceCmd->SetGuidance("Set distance between X and Y views");
  ViewsDistanceCmd->SetParameterName("Size",false);
  ViewsDistanceCmd->SetRange("Size>=0.");
  ViewsDistanceCmd->SetUnitCategory("Length");  
  ViewsDistanceCmd->AvailableForStates(G4State_Idle);

  // calorimeter detector thickness

  CALThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setCALThick",this);
  CALThickCmd->SetGuidance("Set thickness of CAL detectors");
  CALThickCmd->SetParameterName("Size",false);
  CALThickCmd->SetRange("Size>=0.");
  CALThickCmd->SetUnitCategory("Length");  
  CALThickCmd->AvailableForStates(G4State_Idle);

  // number calorimeter detectors 

  NbCALBarsCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALBars",this);
  NbCALBarsCmd->SetGuidance("Set number of CsI Bars.");
  NbCALBarsCmd->SetParameterName("NbSiTiles",false);
  NbCALBarsCmd->SetRange("NbSiTiles>0 && NbSiTiles<100");
  NbCALBarsCmd->AvailableForStates(G4State_Idle);

  // number calorimeter layers

  NbCALLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALLayers",this);
  NbCALLayersCmd->SetGuidance("Set number of CAL Layers.");
  NbCALLayersCmd->SetParameterName("NbCALLayers",false);
  NbCALLayersCmd->SetRange("NbCALLayers>0 && NbCALLayers<16");
  NbCALLayersCmd->AvailableForStates(G4State_Idle);

  // calorimeter detector thickness

  ACDThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setACDThick",this);
  ACDThickCmd->SetGuidance("Set thickness of ACD detectors");
  ACDThickCmd->SetParameterName("Size",false);
  ACDThickCmd->SetRange("Size>=0.");
  ACDThickCmd->SetUnitCategory("Length");  
  ACDThickCmd->AvailableForStates(G4State_Idle);
  
  // update Payload

  UpdateCmd = new G4UIcmdWithoutParameter("/payload/update",this);
  UpdateCmd->SetGuidance("Update payload geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  // magnetic field

  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/payload/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorMessenger::~GammaRayTelDetectorMessenger()
{
  delete ConverterMaterCmd; delete ConverterThickCmd;
  delete NbSiTilesCmd;      delete NbTKRLayersCmd;
  delete SiliconTileXYCmd;  delete SiliconPitchCmd;
  delete SiliconThickCmd;   delete LayerDistanceCmd;
  delete ViewsDistanceCmd;  delete ACDThickCmd;
  delete NbCALLayersCmd;    delete NbCALBarsCmd;
  delete CALThickCmd;       delete UpdateCmd;
  delete MagFieldCmd;       delete GammaRayTeldetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // converter

  if( command == ConverterMaterCmd )
    { GammaRayTelDetector->SetConverterMaterial(newValue);}
   
  if( command == ConverterThickCmd )
    { GammaRayTelDetector->SetConverterThickness(ConverterThickCmd->GetNewDoubleValue(newValue));}
  
  // tracker

  if( command == SiliconTileXYCmd )
    { GammaRayTelDetector->SetTKRTileSizeXY(SiliconTileXYCmd->GetNewDoubleValue(newValue));}

  if( command == SiliconPitchCmd )
    { GammaRayTelDetector->SetTKRSiliconPitch(SiliconPitchCmd->GetNewDoubleValue(newValue));}

  if( command == SiliconThickCmd )
    { GammaRayTelDetector->SetTKRSiliconThickness(SiliconThickCmd->GetNewDoubleValue(newValue));}
  
  if( command == NbSiTilesCmd )
    { GammaRayTelDetector->SetNbOfTKRTiles(NbSiTilesCmd->GetNewIntValue(newValue));}

  if( command == NbTKRLayersCmd )
    { GammaRayTelDetector->SetNbOfTKRLayers(NbTKRLayersCmd->GetNewIntValue(newValue));}

  if( command == LayerDistanceCmd )
    { GammaRayTelDetector->SetTKRLayerDistance(LayerDistanceCmd->GetNewDoubleValue(newValue));}

  if( command == ViewsDistanceCmd )
    { GammaRayTelDetector->SetTKRViewsDistance(ViewsDistanceCmd->GetNewDoubleValue(newValue));}
  
  // calorimeter

  if( command == NbCALLayersCmd )
    { GammaRayTelDetector->SetNbOfCALLayers(NbCALLayersCmd->GetNewIntValue(newValue));}
  
  if( command == NbCALBarsCmd )
    { GammaRayTelDetector->SetNbOfCALBars(NbCALBarsCmd->GetNewIntValue(newValue));}

  if( command == CALThickCmd )
    { GammaRayTelDetector->SetCALBarThickness(CALThickCmd->GetNewDoubleValue(newValue));}
  
  // anticoincidence

  if( command == ACDThickCmd )
    { GammaRayTelDetector->SetACDThickness(ACDThickCmd->GetNewDoubleValue(newValue));}
  
  if( command == UpdateCmd )
    { GammaRayTelDetector->UpdateGeometry(); }
  
  if( command == MagFieldCmd )
    { GammaRayTelDetector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






