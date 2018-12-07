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
 fDetector(Det), fTesthadrDir(0), fDetDir(0), fMaterCmd1(0), fMaterCmd2(0),
 fSizeCmd1(0), fSizeCmd2(0), fSizeCmd3(0),
 fIsotopeCmd(0)
{ 
  fTesthadrDir = new G4UIdirectory("/testhadr/");
  fTesthadrDir->SetGuidance("commands specific to this example");

  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/testhadr/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");

  fMaterCmd1 = new G4UIcmdWithAString("/testhadr/det/setAbsorMat",this);
  fMaterCmd1->SetGuidance("Select absorber material");
  fMaterCmd1->SetParameterName("choice",false);
  fMaterCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMaterCmd2 = new G4UIcmdWithAString("/testhadr/det/setContMat",this);
  fMaterCmd2->SetGuidance("Select container material");
  fMaterCmd2->SetParameterName("choice",false);
  fMaterCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

 fSizeCmd1 = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setAbsorRadius",this);
  fSizeCmd1->SetGuidance("Set absorber radius");
  fSizeCmd1->SetParameterName("Radius",false);
  fSizeCmd1->SetRange("Radius>0.");
  fSizeCmd1->SetUnitCategory("Length");
  fSizeCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

 fSizeCmd2 = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setAbsorLength",this);
  fSizeCmd2->SetGuidance("Set absorber length");
  fSizeCmd2->SetParameterName("Length",false);
  fSizeCmd2->SetRange("Length>0.");
  fSizeCmd2->SetUnitCategory("Length");
  fSizeCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeCmd3 = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setContThick",this);
  fSizeCmd3->SetGuidance("Set container thickness");
  fSizeCmd3->SetParameterName("Thick",false);
  fSizeCmd3->SetRange("Thick>0.");
  fSizeCmd3->SetUnitCategory("Length");
  fSizeCmd3->AvailableForStates(G4State_PreInit,G4State_Idle);
  
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
  delete fMaterCmd1;
  delete fMaterCmd2;
  delete fSizeCmd1;
  delete fSizeCmd2;
  delete fSizeCmd3;
  delete fIsotopeCmd;
  delete fDetDir;
  delete fTesthadrDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd1 )
   { fDetector->SetAbsorMaterial(newValue);}
  if( command == fMaterCmd2 )
   { fDetector->SetContainMaterial(newValue);}

  if( command == fSizeCmd1 )
   { fDetector->SetAbsorRadius(fSizeCmd1->GetNewDoubleValue(newValue));}

  if( command == fSizeCmd2 )
   { fDetector->SetAbsorLength(fSizeCmd2->GetNewDoubleValue(newValue));}

  if( command == fSizeCmd3 )
   { fDetector->SetContainThickness(fSizeCmd3->GetNewDoubleValue(newValue));}

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
