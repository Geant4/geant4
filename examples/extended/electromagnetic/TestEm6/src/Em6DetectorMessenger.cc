// Em6DetectorMessenger.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6DetectorMessenger.hh"

#include "Em6DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6DetectorMessenger::Em6DetectorMessenger(Em6DetectorConstruction * Det)
:Em6Detector(Det)
{
  Em6detDir = new G4UIdirectory("/calor/");
  Em6detDir->SetGuidance("Em6 detector control.");

  MaterCmd = new G4UIcmdWithAString("/calor/setMat",this);
  MaterCmd->SetGuidance("Select Material.");
  MaterCmd->SetParameterName("material",false);
  MaterCmd->AvailableForStates(Idle);

  LBinCmd = new G4UIcmdWith3Vector("/calor/setLbin",this);
  LBinCmd->SetGuidance("set longitudinal bining");
  LBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  LBinCmd->SetParameterName("nLtot","dLradl"," ",true);
  LBinCmd->SetRange("nLtot>=1 && dLradl>0");
  LBinCmd->AvailableForStates(PreInit,Idle);

  FieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  FieldCmd->SetParameterName("Bz",false);
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6DetectorMessenger::~Em6DetectorMessenger()
{
  delete MaterCmd;
  delete LBinCmd;
  delete FieldCmd;
  delete UpdateCmd;
  delete Em6detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == MaterCmd )
   { Em6Detector->SetMaterial(newValue);}

  if( command == LBinCmd )
   { Em6Detector->SetLBining(LBinCmd->GetNew3VectorValue(newValue));}

  if( command == FieldCmd )
   { Em6Detector->SetMagField(FieldCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { Em6Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
