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
// $Id: DetectorMessenger.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(), 
 fDetector(Det), fTestemDir(nullptr), fDetDir(nullptr), fMaterCmd(nullptr),
 fThickCmd(nullptr), fSizeYZCmd(nullptr), fIsotopeCmd(nullptr)
{
  fTestemDir = new G4UIdirectory("/testhadr/");
  fTestemDir->SetGuidance("commands specific to this example");

  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/testhadr/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");

  fMaterCmd = new G4UIcmdWithAString("/testhadr/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fThickCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setThickness",this);
  fThickCmd->SetGuidance("Set thickness of the absor");
  fThickCmd->SetParameterName("Thickness",false);
  fThickCmd->SetRange("Thickness>0.");
  fThickCmd->SetUnitCategory("Length");
  fThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set transverse size of the absor");
  fSizeYZCmd->SetParameterName("Size",false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of density");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm->SetParameterCandidates(unitList);
  fIsotopeCmd->SetParameter(unitPrm);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fThickCmd;
  delete fSizeYZCmd;
  delete fIsotopeCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetAbsorMaterial(newValue);}

  if( command == fThickCmd )
   { fDetector->SetAbsorThickness(fThickCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeYZCmd )
   { fDetector->SetAbsorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}

  if (command == fIsotopeCmd)
   {
     G4int Z; G4int A; G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> Z >> A >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fDetector->MaterialWithSingleIsotope (name,name,dens,Z,A);
     fDetector->SetAbsorMaterial(name);    
   }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
