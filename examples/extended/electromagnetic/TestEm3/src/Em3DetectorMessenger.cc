// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3DetectorMessenger.cc,v 1.2 1999-12-15 14:49:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3DetectorMessenger.hh"

#include "Em3DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3DetectorMessenger::Em3DetectorMessenger(Em3DetectorConstruction * Em3Det)
:Em3Detector(Em3Det)
{ 
  Em3detDir = new G4UIdirectory("/calor/");
  Em3detDir->SetGuidance("Em3 detector control.");
  
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  SizeYZCmd->SetParameterName("Size",false);
  SizeYZCmd->SetRange("Size>0.");
  SizeYZCmd->SetUnitCategory("Length");
  SizeYZCmd->AvailableForStates(Idle);
  
  NbLayersCmd = new G4UIcmdWithAnInteger("/calor/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of layers.");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>0 && NbLayers<500");
  NbLayersCmd->AvailableForStates(Idle);
  
  NbAbsorCmd = new G4UIcmdWithAnInteger("/calor/setNbOfAbsor",this);
  NbAbsorCmd->SetGuidance("Set number of Absorbers.");
  NbAbsorCmd->SetParameterName("NbAbsor",false);
  NbAbsorCmd->SetRange("NbAbsor>0");
  NbAbsorCmd->AvailableForStates(Idle);
   
  AbsorCmd = new G4UIcommand("/calor/setAbsor",this);
  AbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  AbsorCmd->SetGuidance("  absor number : from 0 to NbOfAbsor-1");
  AbsorCmd->SetGuidance("  material name");
  AbsorCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absor number : from 0 to NbOfAbsor-1");
  AbsNbPrm->SetParameterRange("AbsorNb>=0");
  AbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name");
  AbsorCmd->SetParameter(MatPrm);
  //    
  G4UIparameter* ThickPrm = new G4UIparameter("thickness",'d',false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  AbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitCandidates = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitCandidates);
  AbsorCmd->SetParameter(unitPrm);
  //
  AbsorCmd->AvailableForStates(Idle);
  
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(Idle); 
     
  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3DetectorMessenger::~Em3DetectorMessenger()
{
  delete SizeYZCmd;
  delete NbLayersCmd;
  delete NbAbsorCmd;
  delete AbsorCmd;  
  delete MagFieldCmd;
  delete UpdateCmd;
  delete Em3detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    
  if( command == SizeYZCmd )
   { Em3Detector->SetCalorSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == NbLayersCmd )
   { Em3Detector->SetNbOfLayers(NbLayersCmd->GetNewIntValue(newValue));}
   
  if( command == NbAbsorCmd )
   { Em3Detector->SetNbOfAbsor(NbAbsorCmd->GetNewIntValue(newValue));}
   
  if (command == AbsorCmd)
   {
     G4int num; G4double tick;
     char mat[30],unts[30];
     const char* t = newValue;
     G4std::istrstream is((char*)t);
     is >> num >> mat >> tick >> unts;
     G4String material=mat, unt=unts;
     tick *= G4UIcommand::ValueOf(unt);
     Em3Detector->SetAbsorMaterial (num,material);
     Em3Detector->SetAbsorThickness(num,tick);
   }
   
  if( command == MagFieldCmd )
   { Em3Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
           
  if( command == UpdateCmd )
   { Em3Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
