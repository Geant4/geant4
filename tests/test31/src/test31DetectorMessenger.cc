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
// -------------------------------------------------------------
//
// 
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31DetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31DetectorMessenger.hh"
#include "test31DetectorConstruction.hh"
#include "test31Histo.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31DetectorMessenger::test31DetectorMessenger(test31DetectorConstruction* h):
  hDet(h)
{ 
  test31detDir = new G4UIdirectory("/test31/");
  test31detDir->SetGuidance("General test31 commands");
  test31detDir1= new G4UIdirectory("/test31/physics/");
  test31detDir1->SetGuidance("test31 commands to define physics");
  test31detDir2= new G4UIdirectory("/test31/gun/");
  test31detDir2->SetGuidance("test31 commands to define gun");
  if(hDet->GetVerbose() > 0) {
    G4cout << "test31DetectorMessenger: Is constructed" << G4endl;
  }
      
  AbsMaterCmd = new G4UIcmdWithAString("/test31/AbsorberMaterial",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("AbsoberMaterial",false);
  AbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/test31/WorldMaterial",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("WorldMaterial",false);
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/test31/NumberOfAbsorbers",this);
  NumOfAbsCmd->SetGuidance("Set number of absorbers");
  NumOfAbsCmd->SetParameterName("Nabs",false);
  NumOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/test31/AbsorberThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false);  
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->SetUnitCategory("Length");  
  AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbsGapCmd = new G4UIcmdWithADoubleAndUnit("/test31/AbsorberGap",this);
  AbsGapCmd->SetGuidance("Set gap between absorbers");
  AbsGapCmd->SetParameterName("SizeZ",false);  
  AbsGapCmd->SetRange("SizeZ>0.");
  AbsGapCmd->SetUnitCategory("Length");  
  AbsGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/test31/AbsorberXY",this);
  AbsSizYZCmd->SetGuidance("Set sizeXY of the Absorber");
  AbsSizYZCmd->SetParameterName("SizeYZ",false);
  AbsSizYZCmd->SetRange("SizeYZ>0.");
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/test31/WorldZ",this);
  WorldXCmd->SetGuidance("Set Z size of the World");
  WorldXCmd->SetParameterName("WSizeX",false);
  WorldXCmd->SetRange("WSizeX>0.");
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  UpdateCmd = new G4UIcmdWithoutParameter("/test31/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
      
  XMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/test31/FieldX",this);  
  XMagFieldCmd->SetGuidance("Define magnetic field along X");
  XMagFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  XMagFieldCmd->SetParameterName("Bx",false);
  XMagFieldCmd->SetUnitCategory("Magnetic flux density");
  XMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  YMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/test31/FieldY",this);  
  YMagFieldCmd->SetGuidance("Define magnetic field along Y");
  YMagFieldCmd->SetGuidance("Magnetic field will be in Y direction.");
  YMagFieldCmd->SetParameterName("By",false);
  YMagFieldCmd->SetUnitCategory("Magnetic flux density");
  YMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ZMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/test31/FieldZ",this);  
  ZMagFieldCmd->SetGuidance("Define magnetic field along Z");
  ZMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  ZMagFieldCmd->SetParameterName("Bz",false);
  ZMagFieldCmd->SetUnitCategory("Magnetic flux density");
  ZMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  HistoCmd = new G4UIcmdWithAString("/test31/HistoName",this);
  HistoCmd->SetGuidance("Set the name of the histo file");
  HistoCmd->SetParameterName("histo",false);
  HistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ntupCmd = new G4UIcmdWithABool("/hTest/ntuple",this);
  ntupCmd->SetGuidance("Set number ntuple to fill"); 
  ntupCmd->SetParameterName("ntuple",false);
  ntupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfEvt = new G4UIcmdWithAnInteger("/test31/NumberOfEvents",this);
  NumOfEvt->SetGuidance("Set number of event to be simulated");
  NumOfEvt->SetParameterName("Nevt",false);
  NumOfEvt->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/test31/verbose",this);
  verbCmd->SetGuidance("Set verbose for test31");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  intCmd = new G4UIcmdWithAnInteger("/test31/numberAbsToSave",this);
  intCmd->SetGuidance("Set number of absorbers for which "); 
  intCmd->SetGuidance("the energy is saved to tuple");
  intCmd->SetParameterName("numberAbsToSave",false);
  intCmd->AvailableForStates(G4State_PreInit);

  nhistCmd = new G4UIcmdWithAnInteger("/test31/HistoNumber",this);
  nhistCmd->SetGuidance("Set number of histograms to fill"); 
  nhistCmd->SetParameterName("HistoNumber",false);
  nhistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugSCmd = new G4UIcmdWithAnInteger("/test31/nFirstEventToDebug",this);
  nDebugSCmd->SetGuidance("Set number of the first event to debug"); 
  nDebugSCmd->SetParameterName("nFirstEventToDebug",false);
  nDebugSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugECmd = new G4UIcmdWithAnInteger("/test31/nLastEventToDebug",this);
  nDebugECmd->SetGuidance("Set number of the last event to debug"); 
  nDebugECmd->SetParameterName("nLastEventToDebug",false);
  nDebugECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaECmd = new G4UIcmdWithADoubleAndUnit("/test31/maxDeltaEnergy",this);  
  DeltaECmd->SetGuidance("Define scale of delta-Energy histogram");
  DeltaECmd->SetParameterName("DeltaE",false);
  DeltaECmd->SetUnitCategory("Energy");
  DeltaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31DetectorMessenger::~test31DetectorMessenger()
{
  delete NumOfAbsCmd; 
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsGapCmd; 
  delete AbsSizYZCmd;  
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete UpdateCmd;
  delete XMagFieldCmd;
  delete YMagFieldCmd;
  delete ZMagFieldCmd;
  delete HistoCmd;
  delete NumOfEvt;
  delete verbCmd;
  delete intCmd;
  delete nhistCmd;
  delete nDebugSCmd;
  delete nDebugECmd;
  delete test31detDir;
  delete test31detDir1;
  delete test31detDir2;
  delete DeltaECmd;
  delete ntupCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(hDet->GetVerbose() > 1) {
    G4cout << "test31DetectorMessenger: new value = " << newValue << G4endl;
  }

  if( command == NumOfAbsCmd )
   { hDet->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));
     (test31Histo::GetPointer())->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));
   }

  if( command == AbsMaterCmd )
   { hDet->SetAbsorberMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { hDet->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
    (test31Histo::GetPointer())->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
   }

  if( command == AbsGapCmd )
   { hDet->SetGap(AbsGapCmd->GetNewDoubleValue(newValue));
    (test31Histo::GetPointer())->SetGap(AbsGapCmd->GetNewDoubleValue(newValue));
   }

  if( command == WorldMaterCmd )
   { hDet->SetWorldMaterial(newValue);}
   
  if( command == AbsSizYZCmd )
   { hDet->SetAbsorberSizeXY(AbsSizYZCmd->GetNewDoubleValue(newValue));}
      
  if( command == WorldXCmd )
   { hDet->SetWorldSizeZ(WorldXCmd->GetNewDoubleValue(newValue));}
      
  if( command == UpdateCmd )
   { hDet->UpdateGeometry(); }

  if( command == XMagFieldCmd )
   { hDet->SetMagField(XMagFieldCmd->GetNewDoubleValue(newValue),1);}

  if( command == YMagFieldCmd )
   { hDet->SetMagField(YMagFieldCmd->GetNewDoubleValue(newValue),2);}

  if( command == ZMagFieldCmd )
   { hDet->SetMagField(ZMagFieldCmd->GetNewDoubleValue(newValue),3);}

  if( command == HistoCmd ) 
   { (test31Histo::GetPointer())->SetHistoName(newValue);}

  if( command == ntupCmd ) 
   { (test31Histo::GetPointer())->SetNtuple(ntupCmd->GetNewBoolValue(newValue));}

  if( command == NumOfEvt )
   { hDet->SetNumberOfEvents(NumOfAbsCmd->GetNewIntValue(newValue));}

  if( command == verbCmd ){ 
     G4int ver = verbCmd->GetNewIntValue(newValue);
     hDet->SetVerbose(ver);
     (test31Histo::GetPointer())->SetVerbose(ver);
   }

  if( command == intCmd )
   { hDet->SetNumAbsorbersSaved(intCmd->GetNewIntValue(newValue));}

  if( command == nhistCmd )
   { (test31Histo::GetPointer())->SetHistoNumber(nhistCmd->GetNewIntValue(newValue));}

  if( command == nDebugSCmd )
   { hDet->SetFirstEventToDebug(nDebugSCmd->GetNewIntValue(newValue));}

  if( command == nDebugECmd )
   { hDet->SetLastEventToDebug(nDebugECmd->GetNewIntValue(newValue));}

  if( command == DeltaECmd )
   { (test31Histo::GetPointer())
      ->SetMaxEnergy(DeltaECmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
