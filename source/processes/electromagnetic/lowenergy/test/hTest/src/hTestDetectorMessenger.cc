// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
// 
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestDetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestDetectorMessenger.hh"

#include "hTestDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestDetectorMessenger::hTestDetectorMessenger(hTestDetectorConstruction* h):
  hDet(h)
{ 
  hTestdetDir = new G4UIdirectory("/hTest/");
  hTestdetDir->SetGuidance("hTest commands");
      
  AbsMaterCmd = new G4UIcmdWithAString("/hTest/setAbsoberMaterial",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/hTest/setWorldMaterial",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",false);
  WorldMaterCmd->AvailableForStates(Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/hTest/setAbsorberNumber",this);
  NumOfAbsCmd->SetGuidance("Set number of absorbers");
  NumOfAbsCmd->SetParameterName("Nabs",false);
  NumOfAbsCmd->AvailableForStates(Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setAbsorberThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false);  
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->SetUnitCategory("Length");  
  AbsThickCmd->AvailableForStates(Idle);
  
  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setAbsorberXY",this);
  AbsSizYZCmd->SetGuidance("Set sizeXY of the Absorber");
  AbsSizYZCmd->SetParameterName("SizeYZ",false);
  AbsSizYZCmd->SetRange("SizeYZ>0.");
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(Idle);
    
  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setWorldZ",this);
  WorldXCmd->SetGuidance("Set X size of the World");
  WorldXCmd->SetParameterName("WSizeX",false);
  WorldXCmd->SetRange("WSizeX>0.");
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(Idle);
    
  UpdateCmd = new G4UIcmdWithoutParameter("/hTest/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  XMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setFieldX",this);  
  XMagFieldCmd->SetGuidance("Define magnetic field along X");
  XMagFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  XMagFieldCmd->SetParameterName("Bx",false);
  XMagFieldCmd->SetUnitCategory("Magnetic flux density");
  XMagFieldCmd->AvailableForStates(Idle);  

  YMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setFieldY",this);  
  YMagFieldCmd->SetGuidance("Define magnetic field along Y");
  YMagFieldCmd->SetGuidance("Magnetic field will be in Y direction.");
  YMagFieldCmd->SetParameterName("Bx",false);
  YMagFieldCmd->SetUnitCategory("Magnetic flux density");
  YMagFieldCmd->AvailableForStates(Idle);  

  ZMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setFieldZ",this);  
  ZMagFieldCmd->SetGuidance("Define magnetic field along Z");
  ZMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  ZMagFieldCmd->SetParameterName("Bx",false);
  ZMagFieldCmd->SetUnitCategory("Magnetic flux density");
  ZMagFieldCmd->AvailableForStates(Idle);  

  HistoCmd = new G4UIcmdWithAString("/hTest/setHistoName",this);
  HistoCmd->SetGuidance("Set the name of the histo file");
  HistoCmd->SetParameterName("histo",false);
  HistoCmd->AvailableForStates(Idle);

  NumOfEvt = new G4UIcmdWithAnInteger("/hTest/eventNumber",this);
  NumOfEvt->SetGuidance("Set number of event to be simulated");
  NumOfEvt->SetParameterName("Nevt",false);
  NumOfEvt->AvailableForStates(Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestDetectorMessenger::~hTestDetectorMessenger()
{
  delete NumOfAbsCmd; 
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsSizYZCmd;  
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete UpdateCmd;
  delete XMagFieldCmd;
  delete YMagFieldCmd;
  delete ZMagFieldCmd;
  delete HistoCmd;
  delete NumOfEvt;
  delete hTestdetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == NumOfAbsCmd )
   { hDet->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));}

  if( command == AbsMaterCmd )
   { hDet->SetAbsorberMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { hDet->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}

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
   { hDet->SetHistoName(newValue);}

  if( command == NumOfEvt )
   { hDet->SetNumberOfEvents(NumOfAbsCmd->GetNewIntValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
