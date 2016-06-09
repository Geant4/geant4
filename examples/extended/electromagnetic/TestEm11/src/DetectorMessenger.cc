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
// $Id: DetectorMessenger.cc,v 1.2 2006-06-29 16:40:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include <sstream>

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
  
  NbAbsorCmd = new G4UIcmdWithAnInteger("/testem/det/setNbOfAbsor",this);
  NbAbsorCmd->SetGuidance("Set number of Absorbers.");
  NbAbsorCmd->SetParameterName("NbAbsor",false);
  NbAbsorCmd->SetRange("NbAbsor>0");
  NbAbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  AbsorCmd = new G4UIcommand("/testem/det/setAbsor",this);
  AbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  AbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  AbsorCmd->SetGuidance("  material name");
  AbsorCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
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
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  AbsorCmd->SetParameter(unitPrm);
  //
  AbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  NdivCmd = new G4UIcommand("/testem/det/nDivAbsor",this);
  NdivCmd->SetGuidance("Divide the absor nb : number of divisions");
  NdivCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  NdivCmd->SetGuidance("  number of divisions >= 0");
  //
  G4UIparameter* AbsNbPar = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPar->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPar->SetParameterRange("AbsorNb>0");
  NdivCmd->SetParameter(AbsNbPar);
  //
  G4UIparameter* NdivPrm = new G4UIparameter("NdivNb",'i',false);
  NdivPrm->SetGuidance("nb of divisions > 0");
  NdivPrm->SetParameterRange("NdivNb>0");
  NdivCmd->SetParameter(NdivPrm);
  //
  NdivCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
      
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set sizeYZ of the absorber");
  SizeYZCmd->SetParameterName("SizeYZ",false);
  SizeYZCmd->SetRange("SizeYZ>0.");
  SizeYZCmd->SetUnitCategory("Length");
  SizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
            
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
      
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete NbAbsorCmd;
  delete AbsorCmd;
  delete NdivCmd;    
  delete SizeYZCmd;       
  delete MagFieldCmd;
  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    
  if( command == NbAbsorCmd )
   { Detector->SetNbOfAbsor(NbAbsorCmd->GetNewIntValue(newValue));}
   
  if (command == AbsorCmd)
   {
     G4int num; G4double tick;
     G4String unt, mat;
     std::istringstream is(newValue);
     is >> num >> mat >> tick >> unt;
     G4String material=mat;
     tick *= G4UIcommand::ValueOf(unt);
     Detector->SetAbsorMaterial (num,material);
     Detector->SetAbsorThickness(num,tick);
   }
   
  if (command == NdivCmd)
   {
     G4int num, ndiv;
     std::istringstream is(newValue);
     is >> num >> ndiv;
     Detector->SetNbOfDivisions (num,ndiv);
   }
      
  if( command == SizeYZCmd )
   { Detector->SetAbsorSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}   
               
  if( command == MagFieldCmd )
   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
              
  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
