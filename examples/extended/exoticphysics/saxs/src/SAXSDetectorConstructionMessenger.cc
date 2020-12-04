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
/// \file SAXSDetectorConstructionMessenger.cc
/// \brief Implementation of the SAXSDetectorConstructionMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SAXSDetectorConstructionMessenger.hh"
#include "SAXSDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

#include "G4RunManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSDetectorConstructionMessenger::SAXSDetectorConstructionMessenger
(SAXSDetectorConstruction* detconstr):
  G4UImessenger(), fDetector(detconstr)
{
  fCmdDir = new G4UIdirectory("/det/");
  fCmdDir->SetGuidance("Detector Control");

  fSetCustomMatFFfilename = new G4UIcmdWithAString("/det/SetCustomMatFF",this);  
  fSetCustomMatFFfilename->SetGuidance("Set CustomMat FF filename");
  fSetCustomMatFFfilename->SetParameterName("mmff",false);
  fSetCustomMatFFfilename->AvailableForStates(G4State_PreInit);

  fSetCustomMatDensityCmd = new G4UIcmdWithADouble("/det/setCustomMatDensity",
                                                   this);
  fSetCustomMatDensityCmd->SetGuidance("Set density for custom material");
  fSetCustomMatDensityCmd->SetParameterName("cmden",false);
  fSetCustomMatDensityCmd->SetRange("cmden>0.");
  fSetCustomMatDensityCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fSetCustomMatHmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatHmassfract",this);
  fSetCustomMatHmassfractCmd->
    SetGuidance("Set H mass fraction for custom material");
  fSetCustomMatHmassfractCmd->SetParameterName("cmHmf",false);
  fSetCustomMatHmassfractCmd->SetRange("cmHmf>=0.");
  fSetCustomMatHmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatCmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatCmassfract",this);
  fSetCustomMatCmassfractCmd->
    SetGuidance("Set C mass fraction for custom material");
  fSetCustomMatCmassfractCmd->SetParameterName("cmCmf",false);
  fSetCustomMatCmassfractCmd->SetRange("cmCmf>=0.");
  fSetCustomMatCmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatNmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatNmassfract",this);
  fSetCustomMatNmassfractCmd->SetGuidance
    ("Set N mass fraction for custom material");
  fSetCustomMatNmassfractCmd->SetParameterName("cmNmf",false);
  fSetCustomMatNmassfractCmd->SetRange("cmNmf>=0.");
  fSetCustomMatNmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatOmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatOmassfract",this);
  fSetCustomMatOmassfractCmd->SetGuidance
    ("Set O mass fraction for custom material");
  fSetCustomMatOmassfractCmd->SetParameterName("cmOmf",false);
  fSetCustomMatOmassfractCmd->SetRange("cmOmf>=0.");
  fSetCustomMatOmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatNamassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatNamassfract",this);
  fSetCustomMatNamassfractCmd->
    SetGuidance("Set Na mass fraction for custom material");
  fSetCustomMatNamassfractCmd->SetParameterName("cmNamf",false);
  fSetCustomMatNamassfractCmd->SetRange("cmNamf>=0.");
  fSetCustomMatNamassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatPmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatPmassfract",this);
  fSetCustomMatPmassfractCmd->
    SetGuidance("Set P mass fraction for custom material");
  fSetCustomMatPmassfractCmd->SetParameterName("cmPmf",false);
  fSetCustomMatPmassfractCmd->SetRange("cmPmf>=0.");
  fSetCustomMatPmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatSmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatSmassfract",this);
  fSetCustomMatSmassfractCmd->
    SetGuidance("Set S mass fraction for custom material");
  fSetCustomMatSmassfractCmd->SetParameterName("cmSmf",false);
  fSetCustomMatSmassfractCmd->SetRange("cmSmf>=0.");
  fSetCustomMatSmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatClmassfractCmd =
    new G4UIcmdWithADouble("/det/setCustomMatClmassfract",this);
  fSetCustomMatClmassfractCmd->
    SetGuidance("Set Cl mass fraction for custom material");
  fSetCustomMatClmassfractCmd->SetParameterName("cmClmf",false);
  fSetCustomMatClmassfractCmd->SetRange("cmClmf>=0.");
  fSetCustomMatClmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);   

  fSetCustomMatKmassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatKmassfract",this);
  fSetCustomMatKmassfractCmd->
    SetGuidance("Set K mass fraction for custom material");
  fSetCustomMatKmassfractCmd->SetParameterName("cmKmf",false);
  fSetCustomMatKmassfractCmd->SetRange("cmKmf>=0.");
  fSetCustomMatKmassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fSetCustomMatCamassfractCmd = new G4UIcmdWithADouble
    ("/det/setCustomMatCamassfract",this);
  fSetCustomMatCamassfractCmd->
    SetGuidance("Set Ca mass fraction for custom material");
  fSetCustomMatCamassfractCmd->SetParameterName("cmCamf",false);
  fSetCustomMatCamassfractCmd->SetRange("cmCamf>=0.");
  fSetCustomMatCamassfractCmd->AvailableForStates(G4State_PreInit,G4State_Idle);   
    
  fPhantomMaterialCmd = new G4UIcmdWithAnInteger("/det/setPhantomMaterial",this);
  fPhantomMaterialCmd->SetGuidance("Set Phantom material");
  fPhantomMaterialCmd->SetParameterName("PhantomMat",false);
  fPhantomMaterialCmd->SetRange("PhantomMat>=1 && PhantomMat<=30");
  fPhantomMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fPhantomDiameterCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setPhantomDiameter",this);
  fPhantomDiameterCmd->SetGuidance("Set Phantom Diameter");
  fPhantomDiameterCmd->SetParameterName("PhantomDiameter",false);
  fPhantomDiameterCmd->SetUnitCategory("Length");
  fPhantomDiameterCmd->SetRange("PhantomDiameter>0.");
  fPhantomDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
        
  fPhantomHeightCmd = new G4UIcmdWithADoubleAndUnit("/det/setPhantomHeight",this);
  fPhantomHeightCmd->SetGuidance("Set Phantom Thickness");
  fPhantomHeightCmd->SetParameterName("PhantomHeight",false);
  fPhantomHeightCmd->SetUnitCategory("Length");
  fPhantomHeightCmd->SetRange("PhantomHeight>0.");
  fPhantomHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fPhantomZCmd = new G4UIcmdWithADoubleAndUnit("/det/setPhantomZ",this);
  fPhantomZCmd->SetGuidance("Set Phantom Z");
  fPhantomZCmd->SetParameterName("PhantomZ",false);
  fPhantomZCmd->SetUnitCategory("Length");
  fPhantomZCmd->SetRange("PhantomZ>=0.");
  fPhantomZCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
        
  fSetComp0Cmd = new G4UIcmdWithADouble("/det/setComp0",this);
  fSetComp0Cmd->SetGuidance("Set Comp0 for medical material");
  fSetComp0Cmd->SetParameterName("c0",false);
  fSetComp0Cmd->SetRange("c0>=0.");
  fSetComp0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
              
  fSetComp1Cmd = new G4UIcmdWithADouble("/det/setComp1",this);
  fSetComp1Cmd->SetGuidance("Set Comp1 for medical material");
  fSetComp1Cmd->SetParameterName("c1",false);
  fSetComp1Cmd->SetRange("c1>=0.");
  fSetComp1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
              
  fSetComp2Cmd = new G4UIcmdWithADouble("/det/setComp2",this);
  fSetComp2Cmd->SetGuidance("Set Comp2 for medical material");
  fSetComp2Cmd->SetParameterName("c2",false);
  fSetComp2Cmd->SetRange("c2>=0.");
  fSetComp2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSetComp3Cmd = new G4UIcmdWithADouble("/det/setComp3",this);
  fSetComp3Cmd->SetGuidance("Set Comp3 for medical material");
  fSetComp3Cmd->SetParameterName("c3",false);
  fSetComp3Cmd->SetRange("c3>=0.");
  fSetComp3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
    
  fThetaSetupCmd = new G4UIcmdWithADouble("/det/setThetaSetup",this);
  fThetaSetupCmd->SetGuidance("Set theta setup (rad)");
  fThetaSetupCmd->SetParameterName("theta",false);
  fThetaSetupCmd->SetRange("theta>=0.");
  fThetaSetupCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSetSlitsCmd = new G4UIcmdWithABool("/det/setSlits", this);
  fSetSlitsCmd->SetGuidance("set Slits");          
  fSetSlitsCmd->SetParameterName("fIWantSlits",true);
  fSetSlitsCmd->SetDefaultValue(false);  
          
  fSlit1ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/det/setSlit1Thickness",this);
  fSlit1ThicknessCmd->SetGuidance("Set Slit1 Thickness");
  fSlit1ThicknessCmd->SetParameterName("Slit1Thickness",false);
  fSlit1ThicknessCmd->SetUnitCategory("Length");
  fSlit1ThicknessCmd->SetRange("Slit1Thickness>0.");
  fSlit1ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit2ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/det/setSlit2Thickness",this);
  fSlit2ThicknessCmd->SetGuidance("Set Slit2 Thickness");
  fSlit2ThicknessCmd->SetParameterName("Slit2Thickness",false);
  fSlit2ThicknessCmd->SetUnitCategory("Length");
  fSlit2ThicknessCmd->SetRange("Slit2Thickness>0.");
  fSlit2ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit3ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/det/setSlit3Thickness",this);
  fSlit3ThicknessCmd->SetGuidance("Set Slit3 Thickness");
  fSlit3ThicknessCmd->SetParameterName("Slit3Thickness",false);
  fSlit3ThicknessCmd->SetUnitCategory("Length");
  fSlit3ThicknessCmd->SetRange("Slit3Thickness>0.");
  fSlit3ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit4ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/det/setSlit4Thickness",this);
  fSlit4ThicknessCmd->SetGuidance("Set Slit4 Thickness");
  fSlit4ThicknessCmd->SetParameterName("Slit4Thickness",false);
  fSlit4ThicknessCmd->SetUnitCategory("Length");
  fSlit4ThicknessCmd->SetRange("Slit4Thickness>0.");
  fSlit4ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit1DistanceCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit1SampleDistance",this);
  fSlit1DistanceCmd->SetGuidance("Set Slit1-to-Sample Distance");
  fSlit1DistanceCmd->SetParameterName("Slit1dist",false);
  fSlit1DistanceCmd->SetUnitCategory("Length");
  fSlit1DistanceCmd->SetRange("Slit1dist>0.");
  fSlit1DistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit2DistanceCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit2SampleDistance",this);
  fSlit2DistanceCmd->SetGuidance("Set Slit2-to-Sample Distance");
  fSlit2DistanceCmd->SetParameterName("Slit2dist",false);
  fSlit2DistanceCmd->SetUnitCategory("Length");
  fSlit2DistanceCmd->SetRange("Slit2dist>0.");
  fSlit2DistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit3DistanceCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit3SampleDistance",this);
  fSlit3DistanceCmd->SetGuidance("Set Slit3-to-Sample Distance");
  fSlit3DistanceCmd->SetParameterName("Slit3dist",false);
  fSlit3DistanceCmd->SetUnitCategory("Length");
  fSlit3DistanceCmd->SetRange("Slit3dist>0.");
  fSlit3DistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit4DistanceCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit4SampleDistance",this);
  fSlit4DistanceCmd->SetGuidance("Set Slit4-to-Sample Distance");
  fSlit4DistanceCmd->SetParameterName("Slit4dist",false);
  fSlit4DistanceCmd->SetUnitCategory("Length");
  fSlit4DistanceCmd->SetRange("Slit4dist>0.");
  fSlit4DistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit1xApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit1xAperture",this);
  fSlit1xApertureCmd->SetGuidance("Set Slit1 x Aperture");
  fSlit1xApertureCmd->SetParameterName("Slit1xAp",false);
  fSlit1xApertureCmd->SetUnitCategory("Length");
  fSlit1xApertureCmd->SetRange("Slit1xAp>0.");
  fSlit1xApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit2xApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit2xAperture",this);
  fSlit2xApertureCmd->SetGuidance("Set Slit2 x Aperture");
  fSlit2xApertureCmd->SetParameterName("Slit2xAp",false);
  fSlit2xApertureCmd->SetUnitCategory("Length");
  fSlit2xApertureCmd->SetRange("Slit2xAp>0.");
  fSlit2xApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit3xApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit3xAperture",this);
  fSlit3xApertureCmd->SetGuidance("Set Slit3 x Aperture");
  fSlit3xApertureCmd->SetParameterName("Slit3xAp",false);
  fSlit3xApertureCmd->SetUnitCategory("Length");
  fSlit3xApertureCmd->SetRange("Slit3xAp>0.");
  fSlit3xApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit4xApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit4xAperture",this);
  fSlit4xApertureCmd->SetGuidance("Set Slit4 x Aperture");
  fSlit4xApertureCmd->SetParameterName("Slit4xAp",false);
  fSlit4xApertureCmd->SetUnitCategory("Length");
  fSlit4xApertureCmd->SetRange("Slit4xAp>0.");
  fSlit4xApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit1yApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit1yAperture",this);
  fSlit1yApertureCmd->SetGuidance("Set Slit1 y Aperture");
  fSlit1yApertureCmd->SetParameterName("Slit1yAp",false);
  fSlit1yApertureCmd->SetUnitCategory("Length");
  fSlit1yApertureCmd->SetRange("Slit1yAp>0.");
  fSlit1yApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit2yApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit2yAperture",this);
  fSlit2yApertureCmd->SetGuidance("Set Slit2 y Aperture");
  fSlit2yApertureCmd->SetParameterName("Slit2yAp",false);
  fSlit2yApertureCmd->SetUnitCategory("Length");
  fSlit2yApertureCmd->SetRange("Slit2yAp>0.");
  fSlit2yApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit3yApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit3yAperture",this);
  fSlit3yApertureCmd->SetGuidance("Set Slit3 y Aperture");
  fSlit3yApertureCmd->SetParameterName("Slit3yAp",false);
  fSlit3yApertureCmd->SetUnitCategory("Length");
  fSlit3yApertureCmd->SetRange("Slit3yAp>0.");
  fSlit3yApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fSlit4yApertureCmd =
    new G4UIcmdWithADoubleAndUnit("/det/setSlit4yAperture",this);
  fSlit4yApertureCmd->SetGuidance("Set Slit4 y Aperture");
  fSlit4yApertureCmd->SetParameterName("Slit4yAp",false);
  fSlit4yApertureCmd->SetUnitCategory("Length");
  fSlit4yApertureCmd->SetRange("Slit4yAp>0.");
  fSlit4yApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fDetectorThicknessCmd = new G4UIcmdWithADoubleAndUnit
    ("/det/setDetectorThickness",this);
  fDetectorThicknessCmd->SetGuidance("Set Detector Thickness");
  fDetectorThicknessCmd->SetParameterName("DetectorThickness",false);
  fDetectorThicknessCmd->SetUnitCategory("Length");
  fDetectorThicknessCmd->SetRange("DetectorThickness>0.");
  fDetectorThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
       
  fDetectorSizeCmd = new G4UIcmdWithADoubleAndUnit("/det/setDetectorSize",this);
  fDetectorSizeCmd->SetGuidance("Set DetectorSize");
  fDetectorSizeCmd->SetParameterName("scrnsize",false);
  fDetectorSizeCmd->SetUnitCategory("Length");
  fDetectorSizeCmd->SetRange("scrnsize>0.");
  fDetectorSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
  fDetectorDistanceCmd = new G4UIcmdWithADoubleAndUnit
    ("/det/setDetectorSampleDistance",this);
  fDetectorDistanceCmd->SetGuidance("Set Detector Distance");
  fDetectorDistanceCmd->SetParameterName("detDist",false);
  fDetectorDistanceCmd->SetUnitCategory("Length");
  fDetectorDistanceCmd->SetRange("detDist>0.");
  fDetectorDistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSDetectorConstructionMessenger::~SAXSDetectorConstructionMessenger()
{
  delete fCmdDir;
    
  delete fSetCustomMatFFfilename;        
  delete fSetCustomMatDensityCmd;
  delete fSetCustomMatHmassfractCmd;
  delete fSetCustomMatCmassfractCmd;
  delete fSetCustomMatNmassfractCmd;
  delete fSetCustomMatOmassfractCmd;
  delete fSetCustomMatNamassfractCmd;
  delete fSetCustomMatPmassfractCmd;
  delete fSetCustomMatSmassfractCmd;
  delete fSetCustomMatClmassfractCmd;
  delete fSetCustomMatKmassfractCmd;
  delete fSetCustomMatCamassfractCmd;
        
  delete fPhantomMaterialCmd;
  delete fPhantomDiameterCmd;
  delete fPhantomHeightCmd;
  delete fPhantomZCmd; 

  delete fSetComp0Cmd;
  delete fSetComp1Cmd;
  delete fSetComp2Cmd;
  delete fSetComp3Cmd;
    
  delete fThetaSetupCmd;   
    
  delete fSetSlitsCmd;
  delete fSlit1ThicknessCmd;
  delete fSlit2ThicknessCmd;
  delete fSlit3ThicknessCmd;
  delete fSlit4ThicknessCmd;
  delete fSlit1DistanceCmd;
  delete fSlit2DistanceCmd;
  delete fSlit3DistanceCmd;
  delete fSlit4DistanceCmd;
  delete fSlit1xApertureCmd;
  delete fSlit2xApertureCmd;
  delete fSlit3xApertureCmd;
  delete fSlit4xApertureCmd;  
  delete fSlit1yApertureCmd;
  delete fSlit2yApertureCmd;
  delete fSlit3yApertureCmd;
  delete fSlit4yApertureCmd;  

  delete fDetectorThicknessCmd;  
  delete fDetectorSizeCmd;
  delete fDetectorDistanceCmd;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSDetectorConstructionMessenger::SetNewValue(G4UIcommand* command,
                                                    G4String newValue)
{
  if (command == fSetCustomMatFFfilename) 
    fDetector->SetCustomMatFF(newValue);
                
  if (command == fSetCustomMatDensityCmd) 
    fDetector->SetCustomMatDensity
      (fSetCustomMatDensityCmd->GetNewDoubleValue(newValue));
                
  if (command == fSetCustomMatHmassfractCmd) 
    fDetector->SetCustomMatHmassfract
      (fSetCustomMatHmassfractCmd->GetNewDoubleValue(newValue));
        
  if (command == fSetCustomMatCmassfractCmd) 
    fDetector->SetCustomMatCmassfract
      (fSetCustomMatCmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatNmassfractCmd) 
    fDetector->SetCustomMatNmassfract
      (fSetCustomMatNmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatOmassfractCmd) 
    fDetector->SetCustomMatOmassfract
      (fSetCustomMatOmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatNamassfractCmd) 
    fDetector->SetCustomMatNamassfract
      (fSetCustomMatNamassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatPmassfractCmd) 
    fDetector->SetCustomMatPmassfract
      (fSetCustomMatPmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatSmassfractCmd) 
    fDetector->SetCustomMatSmassfract
      (fSetCustomMatSmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatClmassfractCmd) 
    fDetector->SetCustomMatClmassfract
      (fSetCustomMatClmassfractCmd->GetNewDoubleValue(newValue));
  if (command == fSetCustomMatKmassfractCmd) 
    fDetector->SetCustomMatKmassfract
      (fSetCustomMatKmassfractCmd->GetNewDoubleValue(newValue));
  if  (command == fSetCustomMatCamassfractCmd) 
    fDetector->SetCustomMatCamassfract
      (fSetCustomMatCamassfractCmd->GetNewDoubleValue(newValue));
        
  if (command == fPhantomMaterialCmd) {
    fDetector->SetPhantomMaterial
      (fPhantomMaterialCmd->GetNewIntValue(newValue));}
  if (command == fPhantomDiameterCmd) {
    fDetector->SetPhantomDiameter
      (fPhantomDiameterCmd->GetNewDoubleValue(newValue));}
  if (command == fPhantomHeightCmd) {
    fDetector->SetPhantomHeight
      (fPhantomHeightCmd->GetNewDoubleValue(newValue));}
  if (command == fPhantomZCmd) {
    fDetector->SetPhantomZ
      (fPhantomZCmd->GetNewDoubleValue(newValue));}
        
  if (command == fSetComp0Cmd) 
    fDetector->SetComp0(fSetComp0Cmd->GetNewDoubleValue(newValue));
  if (command == fSetComp1Cmd) 
    fDetector->SetComp1(fSetComp1Cmd->GetNewDoubleValue(newValue));
  if (command == fSetComp2Cmd) 
    fDetector->SetComp2(fSetComp2Cmd->GetNewDoubleValue(newValue));
  if (command == fSetComp3Cmd) 
    fDetector->SetComp3(fSetComp3Cmd->GetNewDoubleValue(newValue));
        
  if (command == fThetaSetupCmd)
    fDetector->SetThetaSetup
      (fThetaSetupCmd->GetNewDoubleValue(newValue));          
  if (command == fSetSlitsCmd) 
    fDetector->SetSlits(fSetSlitsCmd->GetNewBoolValue(newValue));        
  if (command == fSlit1ThicknessCmd) 
    fDetector->SetSlit1Thickness
      (fSlit1ThicknessCmd->GetNewDoubleValue(newValue));
  if (command == fSlit2ThicknessCmd)
    fDetector->SetSlit2Thickness
      (fSlit2ThicknessCmd->GetNewDoubleValue(newValue));        
  if (command == fSlit3ThicknessCmd)
    fDetector->SetSlit3Thickness
      (fSlit3ThicknessCmd->GetNewDoubleValue(newValue));
  if (command == fSlit4ThicknessCmd) 
    fDetector->SetSlit4Thickness
      (fSlit4ThicknessCmd->GetNewDoubleValue(newValue));                
  if (command == fSlit1DistanceCmd) 
    fDetector->SetSlit1SampleDistance
      (fSlit1DistanceCmd->GetNewDoubleValue(newValue));
  if (command == fSlit2DistanceCmd) 
    fDetector->SetSlit2SampleDistance
      (fSlit2DistanceCmd->GetNewDoubleValue(newValue));
  if (command == fSlit3DistanceCmd) 
    fDetector->SetSlit3SampleDistance
      (fSlit3DistanceCmd->GetNewDoubleValue(newValue));
  if (command == fSlit4DistanceCmd) 
    fDetector->SetSlit4SampleDistance
      (fSlit4DistanceCmd->GetNewDoubleValue(newValue));
  if (command == fSlit1xApertureCmd) 
    fDetector->SetSlit1xAperture
      (fSlit1xApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit2xApertureCmd) 
    fDetector->SetSlit2xAperture
      (fSlit2xApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit3xApertureCmd) 
    fDetector->SetSlit3xAperture
      (fSlit3xApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit4xApertureCmd) 
    fDetector->SetSlit4xAperture
      (fSlit4xApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit1yApertureCmd) 
    fDetector->SetSlit1yAperture
      (fSlit1yApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit2yApertureCmd) 
    fDetector->SetSlit2yAperture
      (fSlit2yApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit3yApertureCmd) 
    fDetector->SetSlit3yAperture
      (fSlit3yApertureCmd->GetNewDoubleValue(newValue));
  if (command == fSlit4yApertureCmd) 
    fDetector->SetSlit4yAperture
      (fSlit4yApertureCmd->GetNewDoubleValue(newValue));          
  if (command == fDetectorThicknessCmd) 
    fDetector->SetDetectorThickness
      (fDetectorThicknessCmd->GetNewDoubleValue(newValue));
  if (command == fDetectorSizeCmd) 
    fDetector->SetDetectorSize
      (fDetectorSizeCmd->GetNewDoubleValue(newValue));        
  if (command == fDetectorDistanceCmd)
    fDetector->SetDetectorSampleDistance
      (fDetectorDistanceCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

