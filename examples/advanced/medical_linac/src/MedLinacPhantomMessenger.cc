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
// $Id: MedLinacPhantomMessenger.cc,v 1.1 2005/07/03 23:27:37 mpiergen Exp $
//
//  Code developed by: M. Piergentili

#include "MedLinacPhantomMessenger.hh"
#include "MedLinacPhantomSD.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//*********************************************************************

MedLinacPhantomMessenger::MedLinacPhantomMessenger(
                                           MedLinacPhantomSD* MedLinacPhantom)
  :MedLinacP(MedLinacPhantom)
{ 

   
  MedLinacDir = new G4UIdirectory("/PhantomSD/");
  MedLinacDir->SetGuidance("phantom parameters in sensitive detector");

  PhantomDimCmd = new G4UIcmdWithADoubleAndUnit("/PhantomSD/dimension",this);
  PhantomDimCmd->SetGuidance("Set phantom dimension (cm)");
  PhantomDimCmd->SetParameterName("phantomDim",false);
  PhantomDimCmd->SetDefaultUnit( "cm" );
  PhantomDimCmd->SetUnitCategory("Length");
  PhantomDimCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NVoxelsCmd = new G4UIcmdWithAnInteger("/PhantomSD/Nvoxels",this);
  NVoxelsCmd->SetGuidance("Set number of voxels along one axis");
  NVoxelsCmd->SetParameterName("numberOfVoxels",false);
  NVoxelsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


}

//****************************************************************************
MedLinacPhantomMessenger::~MedLinacPhantomMessenger()
{
  delete PhantomDimCmd;
  delete NVoxelsCmd;

  delete MedLinacDir;
  
}

//****************************************************************************

void MedLinacPhantomMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == PhantomDimCmd )
   { MedLinacP->SetPhantomDimension(PhantomDimCmd->GetNewDoubleValue(newValue));}

  if( command == NVoxelsCmd )
   { MedLinacP->SetNumberOfPhantomVoxels(NVoxelsCmd->GetNewIntValue(newValue));  }

}
//****************************************************************************
