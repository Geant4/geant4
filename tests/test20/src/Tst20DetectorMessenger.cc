// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20DetectorMessenger.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20DetectorMessenger ------
//
// ************************************************************

#include "Tst20DetectorMessenger.hh"

#include "Tst20DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20DetectorMessenger::Tst20DetectorMessenger(Tst20DetectorConstruction * Tst20Det)
  :Tst20Detector(Tst20Det)

{ 
  Tst20detDir = new G4UIdirectory("/payload/");
  Tst20detDir->SetGuidance("Tst20 payload control.");
  
  // converter material command
  
  ConverterMaterCmd = new G4UIcmdWithAString("/payload/setConvMat",this);
  ConverterMaterCmd->SetGuidance("Select Material of the Converter.");
  ConverterMaterCmd->SetParameterName("choice",false);
  ConverterMaterCmd->AvailableForStates(Idle);

  // converter thickness command
  
  ConverterThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setConvThick",this);
  ConverterThickCmd->SetGuidance("Set Thickness of the Converter");
  ConverterThickCmd->SetParameterName("Size",false);
  ConverterThickCmd->SetRange("Size>=0.");
  ConverterThickCmd->SetUnitCategory("Length");
  ConverterThickCmd->AvailableForStates(Idle);

  // tracker silicon thickness command

  SiliconThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiThick",this);
  SiliconThickCmd->SetGuidance("Set Thickness of the Silicon");
  SiliconThickCmd->SetParameterName("Size",false);
  SiliconThickCmd->SetRange("Size>=0.");
  SiliconThickCmd->SetUnitCategory("Length");
  SiliconThickCmd->AvailableForStates(Idle);

  // tracker silicon pitch command

  SiliconPitchCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiPitch",this);
  SiliconPitchCmd->SetGuidance("Set Pitch of the Silicon Strips");
  SiliconPitchCmd->SetParameterName("Size",false);
  SiliconPitchCmd->SetRange("Size>=0."); 
  SiliconPitchCmd->SetUnitCategory("Length");
  SiliconPitchCmd->AvailableForStates(Idle);
  
  // tracker silicon tile size command
  
  SiliconTileXYCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setSiTileXY",this);
  SiliconTileXYCmd->SetGuidance("Set XY dimensions of Si Tile");
  SiliconTileXYCmd->SetParameterName("Size",false);
  SiliconTileXYCmd->SetRange("Size>=0.");
  SiliconTileXYCmd->SetUnitCategory("Length");  
  SiliconTileXYCmd->AvailableForStates(Idle);
  
  // tracker number of silicon tiles

  NbSiTilesCmd = new G4UIcmdWithAnInteger("/payload/setNbOfSiTiles",this);
  NbSiTilesCmd->SetGuidance("Set number of Si Tiles.");
  NbSiTilesCmd->SetParameterName("NbSiTiles",false);
  NbSiTilesCmd->SetRange("NbSiTiles>0 && NbSiTiles<200");
  NbSiTilesCmd->AvailableForStates(Idle);

  // tracker number of silicon layers

  NbTKRLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfTKRLayers",this);
  NbTKRLayersCmd->SetGuidance("Set number of TKR Layers.");
  NbTKRLayersCmd->SetParameterName("NbTKRLayers",false);
  NbTKRLayersCmd->SetRange("NbTKRLayers>0 && NbTKRLayers<50");
  NbTKRLayersCmd->AvailableForStates(Idle);

  // tracker layer distance

  LayerDistanceCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setLayerDistance",this);
  LayerDistanceCmd->SetGuidance("Set distance between two layers");
  LayerDistanceCmd->SetParameterName("Size",false);
  LayerDistanceCmd->SetRange("Size>=0.");
  LayerDistanceCmd->SetUnitCategory("Length");  
  LayerDistanceCmd->AvailableForStates(Idle);

  // tracker views distance
  
  TileDistanceCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setViewsDistance",this);
  TileDistanceCmd->SetGuidance("Set distance between X and Y views");
  TileDistanceCmd->SetParameterName("Size",false);
  TileDistanceCmd->SetRange("Size>=0.");
  TileDistanceCmd->SetUnitCategory("Length");  
  TileDistanceCmd->AvailableForStates(Idle);
  /*
  // calorimeter detector thickness

  CALThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setCALThick",this);
  CALThickCmd->SetGuidance("Set thickness of CAL detectors");
  CALThickCmd->SetParameterName("Size",false);
  CALThickCmd->SetRange("Size>=0.");
  CALThickCmd->SetUnitCategory("Length");  
  CALThickCmd->AvailableForStates(Idle);

  // number calorimeter detectors 

  NbCALBarsCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALBars",this);
  NbCALBarsCmd->SetGuidance("Set number of CsI Bars.");
  NbCALBarsCmd->SetParameterName("NbSiTiles",false);
  NbCALBarsCmd->SetRange("NbSiTiles>0 && NbSiTiles<100");
  NbCALBarsCmd->AvailableForStates(Idle);

  // number calorimeter layers

  NbCALLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALLayers",this);
  NbCALLayersCmd->SetGuidance("Set number of CAL Layers.");
  NbCALLayersCmd->SetParameterName("NbCALLayers",false);
  NbCALLayersCmd->SetRange("NbCALLayers>0 && NbCALLayers<16");
  NbCALLayersCmd->AvailableForStates(Idle);
  */

  // calorimeter detector thickness

  ACDThickCmd = new G4UIcmdWithADoubleAndUnit
    ("/payload/setACDThick",this);
  ACDThickCmd->SetGuidance("Set thickness of ACD detectors");
  ACDThickCmd->SetParameterName("Size",false);
  ACDThickCmd->SetRange("Size>=0.");
  ACDThickCmd->SetUnitCategory("Length");  
  ACDThickCmd->AvailableForStates(Idle);
  
  // update Payload

  UpdateCmd = new G4UIcmdWithoutParameter("/payload/update",this);
  UpdateCmd->SetGuidance("Update payload geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  // magnetic field

  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/payload/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20DetectorMessenger::~Tst20DetectorMessenger()
{
  delete ConverterMaterCmd; delete ConverterThickCmd;
  delete NbSiTilesCmd;      delete NbTKRLayersCmd;
  delete SiliconTileXYCmd;  delete SiliconPitchCmd;
  delete SiliconThickCmd;   delete LayerDistanceCmd;
  delete TileDistanceCmd;  delete ACDThickCmd;
  delete UpdateCmd;
  delete MagFieldCmd;       delete Tst20detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // converter

  if( command == ConverterMaterCmd )
    { Tst20Detector->SetConverterMaterial(newValue);}
   
  if( command == ConverterThickCmd )
    { Tst20Detector->SetConverterThickness(ConverterThickCmd->GetNewDoubleValue(newValue));}
  
  // tracker

  if( command == SiliconTileXYCmd )
    { Tst20Detector->SetTKRTileSizeXY(SiliconTileXYCmd->GetNewDoubleValue(newValue));}

  if( command == SiliconPitchCmd )
    { Tst20Detector->SetTKRSiliconPitch(SiliconPitchCmd->GetNewDoubleValue(newValue));}

  if( command == SiliconThickCmd )
    { Tst20Detector->SetTKRSiliconThickness(SiliconThickCmd->GetNewDoubleValue(newValue));}
  
  if( command == NbSiTilesCmd )
    { Tst20Detector->SetNbOfTKRTiles(NbSiTilesCmd->GetNewIntValue(newValue));}

  if( command == NbTKRLayersCmd )
    { Tst20Detector->SetNbOfTKRLayers(NbTKRLayersCmd->GetNewIntValue(newValue));}

  if( command == LayerDistanceCmd )
    { Tst20Detector->SetTKRLayerDistance(LayerDistanceCmd->GetNewDoubleValue(newValue));}

  if( command == TileDistanceCmd )
    { Tst20Detector->SetTKRTileDistance(TileDistanceCmd->GetNewDoubleValue(newValue));}
  
  // anticoincidence
  
  if( command == ACDThickCmd )
    { Tst20Detector->SetACDThickness(ACDThickCmd->GetNewDoubleValue(newValue));}
  
  if( command == ACDMaterialCmd )
    { Tst20Detector->SetACDMaterial(newValue);}
  
  // total

  if( command == UpdateCmd )
    { Tst20Detector->UpdateGeometry(); }
  
  if( command == MagFieldCmd )
    { Tst20Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






