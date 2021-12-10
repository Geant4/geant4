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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
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
:G4UImessenger(),fDetector(Det),
 fTestemDir(nullptr),
 fDetDir(nullptr),
 fSizeYZCmd(nullptr),
 fNbLayersCmd(nullptr),
 fNbAbsorCmd(nullptr),
 fAbsorCmd(nullptr),
 fIsotopeCmd(nullptr)
{ 
  fTestemDir = new G4UIdirectory("/testhadr/");
  fTestemDir->SetGuidance("UI commands specific to this example");
  
  fDetDir = new G4UIdirectory("/testhadr/det/");
  fDetDir->SetGuidance("detector construction commands");
  
  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("Size",false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit);
  fSizeYZCmd->SetToBeBroadcasted(false);
    
  fNbLayersCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfLayers",this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers",false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);
    
  fNbAbsorCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfAbsor",this);
  fNbAbsorCmd->SetGuidance("Set number of Absorbers.");
  fNbAbsorCmd->SetParameterName("NbAbsor",false);
  fNbAbsorCmd->SetRange("NbAbsor>0");
  fNbAbsorCmd->AvailableForStates(G4State_PreInit);
  fNbAbsorCmd->SetToBeBroadcasted(false);
     
  fAbsorCmd = new G4UIcommand("/testhadr/det/setAbsor",this);
  fAbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0."); 
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(MatPrm);
  //    
  G4UIparameter* ThickPrm = new G4UIparameter("thickness",'d',false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  //
  fAbsorCmd->AvailableForStates(G4State_PreInit);
  fAbsorCmd->SetToBeBroadcasted(false);
  
  fIsotopeCmd = new G4UIcommand("/testhadr/det/setIsotopeMat",this);
  fIsotopeCmd->SetGuidance("Build and select a material with single isotope");
  fIsotopeCmd->SetGuidance("  symbol of isotope, Z, A, density of material");
  //
  G4UIparameter* symbPrm = new G4UIparameter("isotope",'s',false);
  symbPrm->SetGuidance("isotope symbol");
  fIsotopeCmd->SetParameter(symbPrm);
  //      
  G4UIparameter* ZPrm = new G4UIparameter("Z",'i',false);
  ZPrm->SetGuidance("Z");
  ZPrm->SetParameterRange("Z>0");
  fIsotopeCmd->SetParameter(ZPrm);
  //      
  G4UIparameter* APrm = new G4UIparameter("A",'i',false);
  APrm->SetGuidance("A");
  APrm->SetParameterRange("A>0");
  fIsotopeCmd->SetParameter(APrm);  
  //    
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fIsotopeCmd->SetParameter(densityPrm);
  //
  G4UIparameter* unitPrm1 = new G4UIparameter("unit",'s',false);
  unitPrm1->SetGuidance("unit of density");
  G4String unitList1 = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm1->SetParameterCandidates(unitList1);
  fIsotopeCmd->SetParameter(unitPrm1);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fSizeYZCmd;
  delete fNbLayersCmd;
  delete fNbAbsorCmd;
  delete fAbsorCmd;
  delete fIsotopeCmd;  
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fSizeYZCmd )
   { fDetector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}

  if( command == fNbLayersCmd )
   { fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));}

  if( command == fNbAbsorCmd )
   { fDetector->SetNbOfAbsor(fNbAbsorCmd->GetNewIntValue(newValue));}
   
  if (command == fAbsorCmd)
   {
     G4int num; G4double tick;
     G4String unt, mat;
     std::istringstream is(newValue);
     is >> num >> mat >> tick >> unt;
     G4String material=mat;
     tick *= G4UIcommand::ValueOf(unt);
     fDetector->SetAbsorMaterial (num,material);
     fDetector->SetAbsorThickness(num,tick);
   }
   
  if (command == fIsotopeCmd)
   {
     G4int Z; G4int A; G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> Z >> A >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fDetector->MaterialWithSingleIsotope (name,name,dens,Z,A);
   }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
