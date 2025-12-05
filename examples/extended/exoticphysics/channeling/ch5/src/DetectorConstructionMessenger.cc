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
// gpaterno, October 2025
//
/// \file DetectorConstructionMessenger.cc
/// \brief Implementation of the DetectorConstruction messenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"

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

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* det):
fDetector(det)
{
    fCmdDir = new G4UIdirectory("/crystal/");
    fCmdDir->SetGuidance("crystal Control");  
    

    fHybridSourceCmd = new G4UIcmdWithABool("/det/isHybridSource", this);
    fHybridSourceCmd->SetGuidance("set if it is a HybridSource");      
    fHybridSourceCmd->SetParameterName("SetHybridSource",true);
    fHybridSourceCmd->SetDefaultValue(false); 
    
    
    fCrystalMaterialCmd = new G4UIcmdWithAString("/crystal/setCrystalMaterial",this);  
    fCrystalMaterialCmd->SetGuidance("Set Crystal Material");
    fCrystalMaterialCmd->SetParameterName("matname",true);
    fCrystalMaterialCmd->SetDefaultValue("Si");
               
    fCrystalSizeCmd = new G4UIcmdWith3VectorAndUnit("/crystal/setCrystalSize",this);
    fCrystalSizeCmd->SetGuidance("Set Crystal size");
    fCrystalSizeCmd->SetParameterName("cryX","cryY","cryZ",false);
    fCrystalSizeCmd->SetUnitCategory("Length");
    fCrystalSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCrystalLatticeCmd = new G4UIcmdWithAString("/crystal/setCrystalLattice",this);  
    fCrystalLatticeCmd->SetGuidance("Set Crystal Lattice");
    fCrystalLatticeCmd->SetParameterName("lattice",false);
    fCrystalLatticeCmd->SetDefaultValue("(111)");    
      
    fCrystalAngleXCmd = new G4UIcmdWithADouble("/crystal/setCrystalAngleX",this);
    fCrystalAngleXCmd->SetGuidance("Set crystal orientation with respet to the beam");
    fCrystalAngleXCmd->SetParameterName("angX",false);
    fCrystalAngleXCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fCrystalAngleYCmd = new G4UIcmdWithADouble("/crystal/setCrystalAngleY",this);
    fCrystalAngleYCmd->SetGuidance("Set crystal orientation with respet to the beam");
    fCrystalAngleYCmd->SetParameterName("angY",false);
    fCrystalAngleYCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fCrystalBendingAngleCmd = 
        new G4UIcmdWithADouble("/crystal/setCrystalBendingAngle",this);
    fCrystalBendingAngleCmd->SetGuidance("Set crystal bending angle");
    fCrystalBendingAngleCmd->SetParameterName("bendingAngle",false);
    fCrystalBendingAngleCmd->SetRange("bendingAngle >= 0");
    fCrystalBendingAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fRadModelCmd = new G4UIcmdWithABool("/crystal/setRadiationModel", this);
    fRadModelCmd->SetGuidance("set Radiation Model");      
    fRadModelCmd->SetParameterName("ActivateRadiationModel",true);
    fRadModelCmd->SetDefaultValue(false);  
    
    fOCeffectsCmd = new G4UIcmdWithABool("/crystal/setOCeffects", this);
    fOCeffectsCmd->SetGuidance("set Oriented Crystal effects");      
    fOCeffectsCmd->SetParameterName("OCeffects",true);
    fOCeffectsCmd->SetDefaultValue(false);
    
    fPotentialPathCmd = new G4UIcmdWithAString("/crystal/setChannelingDataPath",this);
    fPotentialPathCmd->
            SetGuidance("Set the path where to find the available data "
                        "for the G4ChannelingFastSimModel "
                        "if different from G4CHANNELINGDATA");
    fPotentialPathCmd->SetParameterName("channelingDataPath",false);
    fPotentialPathCmd->SetDefaultValue("");
      
      
    fRadiatorConverterSepDistanceCmd = 
        new G4UIcmdWithADoubleAndUnit("/det/setRadiatorConverterSepDistance",this);
    fRadiatorConverterSepDistanceCmd->
        SetGuidance("Set Radiator-Converter Separation Distance");
    fRadiatorConverterSepDistanceCmd->SetParameterName("RCsepDist",false);
    fRadiatorConverterSepDistanceCmd->SetUnitCategory("Length");
    fRadiatorConverterSepDistanceCmd->SetRange("RCsepDist>=0.");
    fRadiatorConverterSepDistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fConverterSizeCmd = new G4UIcmdWith3VectorAndUnit("/det/setConverterSize",this);
    fConverterSizeCmd->SetGuidance("Set Converter size");
    fConverterSizeCmd->SetParameterName("pos3cnvX","pos3cnvY","pos3cnvZ",false);
    fConverterSizeCmd->SetUnitCategory("Length");
    fConverterSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fGranularConverterCmd = new G4UIcmdWithABool("/det/setGranularConverter",this);
    fGranularConverterCmd->SetGuidance("set Collimator");      
    fGranularConverterCmd->SetParameterName("IWantCollimator",true);
    fGranularConverterCmd->SetDefaultValue(false);
      
    fSphereRadiusCmd = new G4UIcmdWithADoubleAndUnit("/det/setSphereRadius",this);
    fSphereRadiusCmd->
        SetGuidance("Set the radius of the spheres composing the Granular Converter");
    fSphereRadiusCmd->SetParameterName("sphr",false);
    fSphereRadiusCmd->SetUnitCategory("Length");
    fSphereRadiusCmd->SetRange("sphr>=0.");
    fSphereRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fConverterMaterialCmd = new G4UIcmdWithAString("/det/setConverterMaterial",this);  
    fConverterMaterialCmd->SetGuidance("Set Converter Material");
    fConverterMaterialCmd->SetParameterName("convmatname",true);
    fConverterMaterialCmd->SetDefaultValue("W"); 
          

    fMagneticFieldCmd = new G4UIcmdWithABool("/det/setMagneticField",this);
    fMagneticFieldCmd->SetGuidance("set Magnetic Field");      
    fMagneticFieldCmd->SetParameterName("IWantMagneticField",true);
    fMagneticFieldCmd->SetDefaultValue(false);  
      
    fFieldValueCmd = new G4UIcmdWithADoubleAndUnit("/det/setMagneticFieldValue",this);
    fFieldValueCmd->SetGuidance("Set Magnetic Field By Value");
    fFieldValueCmd->SetParameterName("By",false);
    fFieldValueCmd->SetUnitCategory("Magnetic flux density");
    fFieldValueCmd->SetRange("By>=0.");
    fFieldValueCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fFieldRegionLengthCmd = 
        new G4UIcmdWithADoubleAndUnit("/det/setMagneticFieldRegionLength",this);
    fFieldRegionLengthCmd->SetGuidance("Set Magnetic Field region length");
    fFieldRegionLengthCmd->SetParameterName("brl",false);
    fFieldRegionLengthCmd->SetUnitCategory("Length");
    fFieldRegionLengthCmd->SetRange("brl>0.");
    fFieldRegionLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      
      
      
    fCollimatorCmd = new G4UIcmdWithABool("/det/setCollimator",this);
    fCollimatorCmd->SetGuidance("set Collimator");      
    fCollimatorCmd->SetParameterName("IWantCollimator",true);
    fCollimatorCmd->SetDefaultValue(false);
            
    fCollimatorApertureCmd = 
        new G4UIcmdWithADoubleAndUnit("/det/setCollimatorAperture",this);
    fCollimatorApertureCmd->SetGuidance("Set Collimator Aperture");
    fCollimatorApertureCmd->SetParameterName("CollAp",false);
    fCollimatorApertureCmd->SetUnitCategory("Length");
    fCollimatorApertureCmd->SetRange("CollAp>0.");
    fCollimatorApertureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCollimatorHoleCmd = new G4UIcmdWithAString("/det/setCollimatorHole",this);  
    fCollimatorHoleCmd->SetGuidance("Set Collimator hole shape");
    fCollimatorHoleCmd->SetParameterName("holeshape",true);
    fCollimatorHoleCmd->SetDefaultValue("squared");
    
    fCollimatorThicknessCmd = 
        new G4UIcmdWithADoubleAndUnit("/det/setCollimatorThickness",this);
    fCollimatorThicknessCmd->SetGuidance("Set Collimator thickenss");
    fCollimatorThicknessCmd->SetParameterName("collth",false);
    fCollimatorThicknessCmd->SetUnitCategory("Length");
    fCollimatorThicknessCmd->SetRange("collth>=0.");
    fCollimatorThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCollimatorSideCmd = new G4UIcmdWithADoubleAndUnit("/det/setCollimatorSide",this);
    fCollimatorSideCmd->SetGuidance("Set Collimator Side size");
    fCollimatorSideCmd->SetParameterName("collside",false);
    fCollimatorSideCmd->SetUnitCategory("Length");
    fCollimatorSideCmd->SetRange("collside>=0.");
    fCollimatorSideCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    fRadiatorCollimatorSepDistanceCmd = 
        new G4UIcmdWithADoubleAndUnit("/det/setRadiatorCollimatorSepDistance",this);
    fRadiatorCollimatorSepDistanceCmd->
        SetGuidance("Set Radiator-Collimator Separation Distance");
    fRadiatorCollimatorSepDistanceCmd->SetParameterName("RCsepDist",false);
    fRadiatorCollimatorSepDistanceCmd->SetUnitCategory("Length");
    fRadiatorCollimatorSepDistanceCmd->SetRange("RCsepDist>=0.");
    fRadiatorCollimatorSepDistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
    
    fVirtualDetectorSizeCmd = 
        new G4UIcmdWith3VectorAndUnit("/crystal/setVirtualDetectorSize",this);
    fVirtualDetectorSizeCmd->SetGuidance("Set VirtualDetector size");
    fVirtualDetectorSizeCmd->SetParameterName("vdX","vdY","vdZ",false);
    fVirtualDetectorSizeCmd->SetUnitCategory("Length");
    fVirtualDetectorSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fScoringCrystalExitCmd = new G4UIcmdWithABool("/det/setScoringCrystalExit",this);
    fScoringCrystalExitCmd->SetGuidance("set IWantScoringCrystalExit");      
    fScoringCrystalExitCmd->SetParameterName("IWantScoringCrystalExit",true);
    fScoringCrystalExitCmd->SetDefaultValue(false);
    fScoringCrystalExitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
    delete fCmdDir;
    
    delete fHybridSourceCmd;
    
    delete fCrystalMaterialCmd;
    delete fCrystalSizeCmd;
    delete fCrystalLatticeCmd;
    delete fCrystalAngleXCmd;
    delete fCrystalAngleYCmd;
    delete fCrystalBendingAngleCmd;
    delete fRadModelCmd;    
    delete fOCeffectsCmd;
    delete fPotentialPathCmd;
    
    delete fRadiatorConverterSepDistanceCmd;
    delete fConverterSizeCmd;
    delete fGranularConverterCmd;
    delete fSphereRadiusCmd;
    delete fConverterMaterialCmd;
    
    delete fMagneticFieldCmd;
    delete fFieldValueCmd;
    delete fFieldRegionLengthCmd;
    
    delete fCollimatorCmd;
    delete fCollimatorApertureCmd;
    delete fCollimatorHoleCmd;
    delete fCollimatorThicknessCmd;
    delete fCollimatorSideCmd;
    delete fRadiatorCollimatorSepDistanceCmd;
    
    delete fVirtualDetectorSizeCmd;
    
    delete fScoringCrystalExitCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fHybridSourceCmd) 
        {fDetector->SetHybridSource(fHybridSourceCmd->GetNewBoolValue(newValue));}    
    
    if (command == fCrystalMaterialCmd) 
        {fDetector->SetCrystalMaterial(newValue);}
    if (command == fCrystalSizeCmd) 
        {fDetector->SetCrystalSize(fCrystalSizeCmd->GetNew3VectorValue(newValue));}    
    if (command == fCrystalLatticeCmd) 
        {fDetector->SetCrystalLattice(newValue);}
    if (command == fCrystalAngleXCmd) 
        {fDetector->SetCrystalAngleX(fCrystalAngleXCmd->GetNewDoubleValue(newValue));}    
    if (command == fCrystalAngleYCmd) 
        {fDetector->SetCrystalAngleY(fCrystalAngleYCmd->GetNewDoubleValue(newValue));}    
    if (command == fCrystalBendingAngleCmd) 
        {fDetector->
            SetCrystalBendingAngle(fCrystalBendingAngleCmd->GetNewDoubleValue(newValue));}        
    if (command == fRadModelCmd) 
        {fDetector->SetRadiationModel(fRadModelCmd->GetNewBoolValue(newValue));}    
    if (command == fOCeffectsCmd) 
        {fDetector->SetOCeffects(fOCeffectsCmd->GetNewBoolValue(newValue));}
    if (command == fPotentialPathCmd) 
        {fDetector->SetPotentialPath(newValue);}  

    if (command == fRadiatorConverterSepDistanceCmd) 
        {fDetector->SetRadiatorConverterSepDistance(
            fRadiatorConverterSepDistanceCmd->GetNewDoubleValue(newValue));}
    if (command == fConverterSizeCmd) 
        {fDetector->SetConverterSize(fConverterSizeCmd->GetNew3VectorValue(newValue));}
    if (command == fGranularConverterCmd) 
        {fDetector->
            SetGranularConverter(fGranularConverterCmd->GetNewBoolValue(newValue));}
    if (command == fSphereRadiusCmd) 
        {fDetector->SetSphereRadius(fSphereRadiusCmd->GetNewDoubleValue(newValue));}    
    if (command == fConverterMaterialCmd) 
        {fDetector->SetConverterMaterial(newValue);}
        
    if (command == fMagneticFieldCmd) 
        {fDetector->SetMagneticField(fMagneticFieldCmd->GetNewBoolValue(newValue));}    
    if (command == fFieldValueCmd) 
        {fDetector->SetFieldValue(fFieldValueCmd->GetNewDoubleValue(newValue));}
    if (command == fFieldRegionLengthCmd) 
        {fDetector->
            SetFieldRegionLength(fFieldRegionLengthCmd->GetNewDoubleValue(newValue));}        
        
    if (command == fCollimatorCmd) 
        {fDetector->SetCollimator(fCollimatorCmd->GetNewBoolValue(newValue));}        
    if (command == fCollimatorApertureCmd) 
        {fDetector->
            SetCollimatorAperture(fCollimatorApertureCmd->GetNewDoubleValue(newValue));}
    if (command == fCollimatorHoleCmd) 
        {fDetector->SetCollimatorHole(newValue);}    
    if (command == fCollimatorThicknessCmd) 
        {fDetector->
            SetCollimatorThickness(fCollimatorThicknessCmd->GetNewDoubleValue(newValue));}
    if (command == fCollimatorSideCmd) 
        {fDetector->SetCollimatorSide(fCollimatorSideCmd->GetNewDoubleValue(newValue));}
    if (command == fRadiatorCollimatorSepDistanceCmd) 
        {fDetector->SetRadiatorCollimatorSepDistance(
            fRadiatorCollimatorSepDistanceCmd->GetNewDoubleValue(newValue));}

    if (command == fVirtualDetectorSizeCmd) 
        {fDetector->
            SetVirtualDetectorSize(fVirtualDetectorSizeCmd->GetNew3VectorValue(newValue));}
    if (command == fScoringCrystalExitCmd) 
        {fDetector->
            SetScoringCrystalExit(fScoringCrystalExitCmd->GetNewBoolValue(newValue));}             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

