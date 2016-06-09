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
// $Id: MedLinacPhantomMessenger.cc,v 1.2 2006/06/29 16:04:33 gunter Exp $
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
