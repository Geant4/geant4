//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst51DetectorMessenger.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// 
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "Tst51DetectorMessenger.hh"
#include "Tst51DetectorConstruction.hh"

Tst51DetectorMessenger::Tst51DetectorMessenger(
                                           Tst51DetectorConstruction* tst51Det)
:detector(tst51Det)
{ 
  tst51Dir = new G4UIdirectory("/target/");
  tst51Dir -> SetGuidance("Tst51 target control.");

  targetMaterialCmd = new G4UIcmdWithAString("/target/setMaterial",this);
  targetMaterialCmd -> SetGuidance("Select Material of the Target.");
  targetMaterialCmd -> SetParameterName("choice",false);
  targetMaterialCmd -> AvailableForStates(G4State_Idle);
     
  targetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/target/setThickness",this);
  targetThicknessCmd -> SetGuidance("Set Thickness of the target");
  targetThicknessCmd -> SetParameterName("Size",false);
  targetThicknessCmd -> SetRange("Size>=0.");
  targetThicknessCmd -> SetUnitCategory("Length");
  targetThicknessCmd -> AvailableForStates(G4State_Idle);
       
 }

Tst51DetectorMessenger::~Tst51DetectorMessenger()
{
  delete targetThicknessCmd; 
  delete targetMaterialCmd; 
  delete tst51Dir;
}

void Tst51DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(command == targetMaterialCmd) detector -> SetTargetMaterial(newValue);
   
  if(command == targetThicknessCmd)
  detector -> SetTargetThickness(targetThicknessCmd
                                               ->GetNewDoubleValue(newValue));
  
}
