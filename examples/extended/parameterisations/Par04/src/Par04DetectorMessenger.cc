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
#include "Par04DetectorMessenger.hh"
#include <CLHEP/Units/SystemOfUnits.h>   // for pi
#include <G4ApplicationState.hh>         // for G4State_PreInit, G4State_Idle
#include <G4ThreeVector.hh>              // for G4ThreeVector
#include <G4Types.hh>                    // for G4bool, G4double, G4int
#include <G4UIcommand.hh>                // for G4UIcommand
#include <G4UImessenger.hh>              // for G4UImessenger
#include <G4UIparameter.hh>              // for G4UIparameter
#include <istream>                       // for basic_istream, basic_istream...
#include <string>                        // for operator>>
#include "G4UIcmdWithADoubleAndUnit.hh"  // for G4UIcmdWithADoubleAndUnit
#include "G4UIcmdWithAnInteger.hh"       // for G4UIcmdWithAnInteger
#include "G4UIcmdWithoutParameter.hh"    // for G4UIcmdWithoutParameter
#include "G4UIdirectory.hh"              // for G4UIdirectory
#include "Par04DetectorConstruction.hh"  // for Par04DetectorConstruction

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorMessenger::Par04DetectorMessenger(Par04DetectorConstruction* aDetector)
  : G4UImessenger()
  , fDetector(aDetector)
{
  fExampleDir = new G4UIdirectory("/Par04/");
  fExampleDir->SetGuidance("UI commands specific to this example");

  fDetectorDir = new G4UIdirectory("/Par04/detector/");
  fDetectorDir->SetGuidance("Detector construction UI commands");

  fPrintCmd = new G4UIcmdWithoutParameter("/Par04/detector/print", this);
  fPrintCmd->SetGuidance("Print current settings.");

  fDetectorInnerRadiusCmd =
    new G4UIcmdWithADoubleAndUnit("/Par04/detector/setDetectorInnerRadius", this);
  fDetectorInnerRadiusCmd->SetGuidance("Set cylindrical detector inner radius");
  fDetectorInnerRadiusCmd->SetParameterName("Size", false);
  fDetectorInnerRadiusCmd->SetRange("Size>0.");
  fDetectorInnerRadiusCmd->SetUnitCategory("Length");
  fDetectorInnerRadiusCmd->AvailableForStates(G4State_PreInit);
  fDetectorInnerRadiusCmd->SetToBeBroadcasted(false);

  fDetectorLengthCmd = new G4UIcmdWithADoubleAndUnit("/Par04/detector/setDetectorLength", this);
  fDetectorLengthCmd->SetGuidance("Set length of the detector (cylinder length)");
  fDetectorLengthCmd->SetParameterName("Size", false);
  fDetectorLengthCmd->SetRange("Size>0.");
  fDetectorLengthCmd->SetUnitCategory("Length");
  fDetectorLengthCmd->AvailableForStates(G4State_PreInit);
  fDetectorLengthCmd->SetToBeBroadcasted(false);

  fNbLayersCmd = new G4UIcmdWithAnInteger("/Par04/detector/setNbOfLayers", this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers", false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);

  fAbsorCmd = new G4UIcommand("/Par04/detector/setAbsorber", this);
  fAbsorCmd->SetGuidance("Set the absorber id, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absorber number : from 0 to 1");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0");
  fAbsorCmd->SetGuidance("  if sensitive : true/false.");
  auto  absNbPrm = new G4UIparameter("AbsorNb", 'i', false);
  absNbPrm->SetGuidance("absor number : from 0 to 1");
  absNbPrm->SetParameterRange("AbsorNb>-1&AbsoNb<2");
  fAbsorCmd->SetParameter(absNbPrm);
  auto  matPrm = new G4UIparameter("material", 's', false);
  matPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(matPrm);
  auto  thickPrm = new G4UIparameter("thickness", 'd', false);
  thickPrm->SetGuidance("thickness of absorber");
  thickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(thickPrm);
  auto  unitPrm = new G4UIparameter("unit", 's', false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  auto  sensitivePrm = new G4UIparameter("sensitive", 'b', false);
  sensitivePrm->SetGuidance("if absorber is sensitive (registers energy deposits)");
  fAbsorCmd->SetParameter(sensitivePrm);

  fAbsorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fAbsorCmd->SetToBeBroadcasted(false);

  fMeshDir = new G4UIdirectory("/Par04/mesh/");
  fMeshDir->SetGuidance("Mesh UI commands");

  fMeshNbRhoCellsCmd = new G4UIcmdWithAnInteger("/Par04/mesh/setNbOfRhoCells", this);
  fMeshNbRhoCellsCmd->SetGuidance("Set number of rho cells in the cylindrical mesh readout.");
  fMeshNbRhoCellsCmd->SetParameterName("NbRhoCells", false);
  fMeshNbRhoCellsCmd->SetRange("NbRhoCells>0");
  fMeshNbRhoCellsCmd->AvailableForStates(G4State_PreInit);
  fMeshNbRhoCellsCmd->SetToBeBroadcasted(false);

  fMeshNbPhiCellsCmd = new G4UIcmdWithAnInteger("/Par04/mesh/setNbOfPhiCells", this);
  fMeshNbPhiCellsCmd->SetGuidance("Set number of phi cells in the cylindrical mesh readout.");
  fMeshNbPhiCellsCmd->SetParameterName("NbPhiCells", false);
  fMeshNbPhiCellsCmd->SetRange("NbPhiCells>0");
  fMeshNbPhiCellsCmd->AvailableForStates(G4State_PreInit);
  fMeshNbPhiCellsCmd->SetToBeBroadcasted(false);

  fMeshNbZCellsCmd = new G4UIcmdWithAnInteger("/Par04/mesh/setNbOfZCells", this);
  fMeshNbZCellsCmd->SetGuidance("Set number of z cells in the cylindrical mesh readout.");
  fMeshNbZCellsCmd->SetParameterName("NbZCells", false);
  fMeshNbZCellsCmd->SetRange("NbZCells>0");
  fMeshNbZCellsCmd->AvailableForStates(G4State_PreInit);
  fMeshNbZCellsCmd->SetToBeBroadcasted(false);

  fMeshSizeRhoCellsCmd = new G4UIcmdWithADoubleAndUnit("/Par04/mesh/setSizeOfRhoCells", this);
  fMeshSizeRhoCellsCmd->SetGuidance("Set size of rho cells in the cylindrical readout mesh");
  fMeshSizeRhoCellsCmd->SetParameterName("Size", false);
  fMeshSizeRhoCellsCmd->SetRange("Size>0.");
  fMeshSizeRhoCellsCmd->SetUnitCategory("Length");
  fMeshSizeRhoCellsCmd->AvailableForStates(G4State_PreInit);
  fMeshSizeRhoCellsCmd->SetToBeBroadcasted(false);

  fMeshSizeZCellsCmd = new G4UIcmdWithADoubleAndUnit("/Par04/mesh/setSizeOfZCells", this);
  fMeshSizeZCellsCmd->SetGuidance("Set size of z cells in the cylindrical readout mesh");
  fMeshSizeZCellsCmd->SetParameterName("Size", false);
  fMeshSizeZCellsCmd->SetRange("Size>0.");
  fMeshSizeZCellsCmd->SetUnitCategory("Length");
  fMeshSizeZCellsCmd->AvailableForStates(G4State_PreInit);
  fMeshSizeZCellsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorMessenger::~Par04DetectorMessenger()
{
  delete fPrintCmd;
  delete fDetectorInnerRadiusCmd;
  delete fDetectorLengthCmd;
  delete fNbLayersCmd;
  delete fAbsorCmd;
  delete fDetectorDir;
  delete fMeshNbRhoCellsCmd;
  delete fMeshNbPhiCellsCmd;
  delete fMeshNbZCellsCmd;
  delete fMeshSizeRhoCellsCmd;
  delete fMeshSizeZCellsCmd;
  delete fMeshDir;
  delete fExampleDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorMessenger::SetNewValue(G4UIcommand* aCommand, G4String aNewValue)
{
  if(aCommand == fPrintCmd)
  {
    fDetector->Print();
  }
  else if(aCommand == fDetectorInnerRadiusCmd)
  {
    fDetector->SetInnerRadius(fDetectorInnerRadiusCmd->GetNewDoubleValue(aNewValue));
  }
  else if(aCommand == fDetectorLengthCmd)
  {
    fDetector->SetLength(fDetectorInnerRadiusCmd->GetNewDoubleValue(aNewValue));
  }
  else if(aCommand == fNbLayersCmd)
  {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fAbsorCmd)
  {
    G4int num;
    G4double thick;
    G4String unt, mat;
    G4bool sensitive;
    std::istringstream is(aNewValue);
    is >> num >> mat >> thick >> unt >> std::boolalpha >> sensitive;
    G4String material = mat;
    thick *= G4UIcommand::ValueOf(unt);
    fDetector->SetAbsorberMaterial(num, material);
    fDetector->SetAbsorberThickness(num, thick);
    fDetector->SetAbsorberSensitivity(num, sensitive);
  }
  else if(aCommand == fMeshNbRhoCellsCmd)
  {
    fDetector->SetMeshNbOfCells(0, fMeshNbRhoCellsCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fMeshNbPhiCellsCmd)
  {
    fDetector->SetMeshNbOfCells(1, fMeshNbPhiCellsCmd->GetNewIntValue(aNewValue));
    fDetector->SetMeshSizeOfCells(1,
                                  2. * CLHEP::pi / fMeshNbPhiCellsCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fMeshNbZCellsCmd)
  {
    fDetector->SetMeshNbOfCells(2, fMeshNbZCellsCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fMeshSizeRhoCellsCmd)
  {
    fDetector->SetMeshSizeOfCells(0, fMeshSizeRhoCellsCmd->GetNewDoubleValue(aNewValue));
  }
  else if(aCommand == fMeshSizeZCellsCmd)
  {
    fDetector->SetMeshSizeOfCells(2, fMeshSizeZCellsCmd->GetNewDoubleValue(aNewValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par04DetectorMessenger::GetCurrentValue(G4UIcommand* aCommand)
{
  G4String cv;

  if(aCommand == fDetectorInnerRadiusCmd)
  {
    cv = fDetectorInnerRadiusCmd->ConvertToString(fDetector->GetInnerRadius(), "mm");
  }
  else if(aCommand == fDetectorLengthCmd)
  {
    cv = fDetectorLengthCmd->ConvertToString(fDetector->GetLength(), "mm");
  }
  else if(aCommand == fNbLayersCmd)
  {
    cv = fNbLayersCmd->ConvertToString(fDetector->GetNbOfLayers());
  }
  else if(aCommand == fMeshNbRhoCellsCmd)
  {
    cv = fMeshNbRhoCellsCmd->ConvertToString(fDetector->GetMeshNbOfCells()[0]);
  }
  else if(aCommand == fMeshNbPhiCellsCmd)
  {
    cv = fMeshNbPhiCellsCmd->ConvertToString(fDetector->GetMeshNbOfCells()[1]);
  }
  else if(aCommand == fMeshNbZCellsCmd)
  {
    cv = fMeshNbZCellsCmd->ConvertToString(fDetector->GetMeshNbOfCells()[2]);
  }
  else if(aCommand == fMeshSizeRhoCellsCmd)
  {
    cv = fMeshSizeRhoCellsCmd->ConvertToString(fDetector->GetMeshSizeOfCells()[0]);
  }
  else if(aCommand == fMeshSizeZCellsCmd)
  {
    cv = fMeshSizeZCellsCmd->ConvertToString(fDetector->GetMeshSizeOfCells()[2]);
  }
  return cv;
}
