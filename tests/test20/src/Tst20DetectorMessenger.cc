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
//
// $Id: Tst20DetectorMessenger.cc,v 1.6 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $ 


#include "Tst20DetectorMessenger.hh"

#include "Tst20DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


Tst20DetectorMessenger::Tst20DetectorMessenger(Tst20DetectorConstruction* detConstr) : detector(detConstr)
{ 
  directory = new G4UIdirectory("/calor/");
  directory->SetGuidance("Tst20 detector control.");
      
  absMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  absMaterCmd->SetGuidance("Select Material of the Absorber.");
  absMaterCmd->SetParameterName("choice",true);
  absMaterCmd->SetDefaultValue("Lead");
  absMaterCmd->AvailableForStates(G4State_Idle);
  
  worldMaterCmd = new G4UIcmdWithAString("/calor/setWorldMat",this);
  worldMaterCmd->SetGuidance("Select Material of the World.");
  worldMaterCmd->SetParameterName("wchoice",true);
  worldMaterCmd->SetDefaultValue("Air");
  worldMaterCmd->AvailableForStates(G4State_Idle);
  
  absThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsThick",this);
  absThickCmd->SetGuidance("Set Thickness of the Absorber");
  absThickCmd->SetParameterName("SizeZ",false,false);
  absThickCmd->SetDefaultUnit("mm");
  absThickCmd->SetRange("SizeZ>0.");
  absThickCmd->AvailableForStates(G4State_Idle);
  
  absRadCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsRad",this);
  absRadCmd->SetGuidance("Set radius of the Absorber");
  absRadCmd->SetParameterName("SizeR",false,false);
  absRadCmd->SetDefaultUnit("mm");
  absRadCmd->SetRange("SizeR>0.");
  absRadCmd->AvailableForStates(G4State_Idle);
  
  absZposCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsZpos",this);
  absZposCmd->SetGuidance("Set Z position of the Absorber");
  absZposCmd->SetParameterName("Zpos",false,false);
  absZposCmd->SetDefaultUnit("mm");
  absZposCmd->AvailableForStates(G4State_Idle);
  
  worldZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldZ",this);
  worldZCmd->SetGuidance("Set Z size of the World");
  worldZCmd->SetParameterName("WSizeZ",false,false);
  worldZCmd->SetDefaultUnit("mm");
  worldZCmd->SetRange("WSizeZ>0.");
  worldZCmd->AvailableForStates(G4State_Idle);
  
  worldRCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldR",this);
  worldRCmd->SetGuidance("Set R size of the World");
  worldRCmd->SetParameterName("WSizeR",false,false);
  worldRCmd->SetDefaultUnit("mm");
  worldRCmd->SetRange("WSizeR>0.");
  worldRCmd->AvailableForStates(G4State_Idle);
  
  updateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  updateCmd->SetGuidance("Update calorimeter geometry.");
  updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateCmd->SetGuidance("if you changed geometrical value(s)");
  updateCmd->AvailableForStates(G4State_Idle);
      
  magFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  magFieldCmd->SetGuidance("Define magnetic field");
  magFieldCmd->SetGuidance("Magnetic field will be in Z direction");
  magFieldCmd->SetParameterName("Bz",false,false);
  magFieldCmd->SetDefaultUnit("tesla");
  magFieldCmd->AvailableForStates(G4State_Idle);  

}



Tst20DetectorMessenger::~Tst20DetectorMessenger()
{
  delete absMaterCmd; 
  delete absThickCmd; 
  delete absRadCmd;  
  delete absZposCmd; 
  delete worldMaterCmd;
  delete worldZCmd;
  delete worldRCmd;
  delete updateCmd;
  delete magFieldCmd;
  delete directory;
}



void Tst20DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == absMaterCmd )
    { 
      detector->SetAbsorberMaterial(newValue);
    }
   
  if ( command == worldMaterCmd )
    { 
      detector->SetWorldMaterial(newValue);
    }
   
  if ( command == absThickCmd )
    { 
      detector->SetAbsorberThickness(absThickCmd->GetNewDoubleValue(newValue));
    }
   
  if ( command == absRadCmd )
    { 
      detector->SetAbsorberRadius(absRadCmd->GetNewDoubleValue(newValue));
    }
   
  if ( command == absZposCmd )
    { 
      detector->SetAbsorberZpos(absZposCmd->GetNewDoubleValue(newValue));
    }
   
  if ( command == worldZCmd )
    { 
      detector->SetWorldSizeZ(worldZCmd->GetNewDoubleValue(newValue));
    }
   
  if ( command == worldRCmd )
    { 
      detector->SetWorldSizeR(worldRCmd->GetNewDoubleValue(newValue));
    }
   
  if ( command == updateCmd )
    { 
      detector->UpdateGeometry(); 
    }

  if ( command == magFieldCmd )
    { 
      detector->SetMagField(magFieldCmd->GetNewDoubleValue(newValue));
    }
}


