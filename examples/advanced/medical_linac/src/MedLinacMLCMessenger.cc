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
// $Id: MedLinacMLCMessenger.cc,v 1.1 2004-11-24 16:53:29 mpiergen Exp $
//
//  Code developed by: M. Piergentili

#include "MedLinacMLCMessenger.hh"

#include "MedLinacMLCDecorator.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//*********************************************************************

MedLinacMLCMessenger::MedLinacMLCMessenger(
                                           MedLinacMLCDecorator* MedLinacMLC)
:MedLinacMLCDeco(MedLinacMLC)
{ 

  MedLinacMLCDir = new G4UIdirectory("/MLC/");
  MedLinacMLCDir->SetGuidance("A1 position");
  
  aDir = new G4UIdirectory("/MLC/A/");
  
  a1yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A1position",this);
  a1yPosCmd->SetGuidance("Set A1 position ");
  a1yPosCmd->SetParameterName("A1Pos_y",false);
  a1yPosCmd->SetRange("A1Pos_y>=0. && A1Pos_y<15.0");
  a1yPosCmd->SetDefaultUnit( "cm" );
  a1yPosCmd->SetUnitCategory("Length");
  a1yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a2yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A2position",this);
  a2yPosCmd->SetGuidance("Set A2 position ");
  a2yPosCmd->SetParameterName("A2Pos_y",false);
  a2yPosCmd->SetRange("A2Pos_y>=0. && A2Pos_y<15.0");
  a2yPosCmd->SetDefaultUnit( "cm" );
  a2yPosCmd->SetUnitCategory("Length");
  a2yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  a3yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A3position",this);
  a3yPosCmd->SetGuidance("Set A3 position ");
  a3yPosCmd->SetParameterName("A3Pos_y",false);
  a3yPosCmd->SetRange("A3Pos_y>=0. && A3Pos_y<15.0");
  a3yPosCmd->SetDefaultUnit( "cm" );
  a3yPosCmd->SetUnitCategory("Length");
  a3yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a4yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A4position",this);
  a4yPosCmd->SetGuidance("Set A4 position ");
  a4yPosCmd->SetParameterName("A4Pos_y",false);
  a4yPosCmd->SetRange("A4Pos_y>=0. && A4Pos_y<15.0");
  a4yPosCmd->SetDefaultUnit( "cm" );
  a4yPosCmd->SetUnitCategory("Length");
  a4yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a5yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A5position",this);
  a5yPosCmd->SetGuidance("Set A5 position ");
  a5yPosCmd->SetParameterName("A5Pos_y",false);
  a5yPosCmd->SetRange("A5Pos_y>=0. && A5Pos_y<15.0");
  a5yPosCmd->SetDefaultUnit( "cm" );
  a5yPosCmd->SetUnitCategory("Length");
  a5yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a6yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A6position",this);
  a6yPosCmd->SetGuidance("Set A6 position ");
  a6yPosCmd->SetParameterName("A6Pos_y",false);
  a6yPosCmd->SetRange("A6Pos_y>=0. && A6Pos_y<15.0");
  a6yPosCmd->SetDefaultUnit( "cm" );
  a6yPosCmd->SetUnitCategory("Length");
  a6yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a7yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A7position",this);
  a7yPosCmd->SetGuidance("Set A7 position ");
  a7yPosCmd->SetParameterName("A7Pos_y",false);
  a7yPosCmd->SetRange("A7Pos_y>=0. && A7Pos_y<15.0");
  a7yPosCmd->SetDefaultUnit( "cm" );
  a7yPosCmd->SetUnitCategory("Length");
  a7yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a8yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A8position",this);
  a8yPosCmd->SetGuidance("Set A8 position ");
  a8yPosCmd->SetParameterName("A8Pos_y",false);
  a8yPosCmd->SetRange("A8Pos_y>=0. && A8Pos_y<15.0");
  a8yPosCmd->SetDefaultUnit( "cm" );
  a8yPosCmd->SetUnitCategory("Length");
  a8yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a9yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A9position",this);
  a9yPosCmd->SetGuidance("Set A9 position ");
  a9yPosCmd->SetParameterName("A9Pos_y",false);
  a9yPosCmd->SetRange("A9Pos_y>=0. && A9Pos_y<15.0");
  a9yPosCmd->SetDefaultUnit( "cm" );
  a9yPosCmd->SetUnitCategory("Length");
  a9yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a10yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A10position",this);
  a10yPosCmd->SetGuidance("Set A10 position ");
  a10yPosCmd->SetParameterName("A10Pos_y",false);
  a10yPosCmd->SetRange("A10Pos_y>=0. && A10Pos_y<15.0");
  a10yPosCmd->SetDefaultUnit( "cm" );
  a10yPosCmd->SetUnitCategory("Length");
  a10yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a11yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A11position",this);
  a11yPosCmd->SetGuidance("Set A11 position ");
  a11yPosCmd->SetParameterName("A11Pos_y",false);
  a11yPosCmd->SetRange("A11Pos_y>=0. && A11Pos_y<15.0");
  a11yPosCmd->SetDefaultUnit( "cm" );
  a11yPosCmd->SetUnitCategory("Length");
  a11yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a12yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A12position",this);
  a12yPosCmd->SetGuidance("Set A12 position ");
  a12yPosCmd->SetParameterName("A12Pos_y",false);
  a12yPosCmd->SetRange("A12Pos_y>=0. && A12Pos_y<15.0");
  a12yPosCmd->SetDefaultUnit( "cm" );
  a12yPosCmd->SetUnitCategory("Length");
  a12yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a13yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A13position",this);
  a13yPosCmd->SetGuidance("Set A13 position ");
  a13yPosCmd->SetParameterName("A13Pos_y",false);
  a13yPosCmd->SetRange("A13Pos_y>=0. && A13Pos_y<15.0");
  a13yPosCmd->SetDefaultUnit( "cm" );
  a13yPosCmd->SetUnitCategory("Length");
  a13yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a14yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A14position",this);
  a14yPosCmd->SetGuidance("Set A14 position ");
  a14yPosCmd->SetParameterName("A14Pos_y",false);
  a14yPosCmd->SetRange("A14Pos_y>=0. && A14Pos_y<15.0");
  a14yPosCmd->SetDefaultUnit( "cm" );
  a14yPosCmd->SetUnitCategory("Length");
  a14yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a15yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A15position",this);
  a15yPosCmd->SetGuidance("Set A15 position ");
  a15yPosCmd->SetParameterName("A15Pos_y",false);
  a15yPosCmd->SetRange("A15Pos_y>=0. && A15Pos_y<15.0");
  a15yPosCmd->SetDefaultUnit( "cm" );
  a15yPosCmd->SetUnitCategory("Length");
  a15yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a16yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A16position",this);
  a16yPosCmd->SetGuidance("Set A16 position ");
  a16yPosCmd->SetParameterName("A16Pos_y",false);
  a16yPosCmd->SetRange("A16Pos_y>=0. && A16Pos_y<15.0");
  a16yPosCmd->SetDefaultUnit( "cm" );
  a16yPosCmd->SetUnitCategory("Length");
  a16yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a17yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A17position",this);
  a17yPosCmd->SetGuidance("Set A17 position ");
  a17yPosCmd->SetParameterName("A17Pos_y",false);
  a17yPosCmd->SetRange("A17Pos_y>=0. && A17Pos_y<15.0");
  a17yPosCmd->SetDefaultUnit( "cm" );
  a17yPosCmd->SetUnitCategory("Length");
  a17yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a18yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A18position",this);
  a18yPosCmd->SetGuidance("Set A18 position ");
  a18yPosCmd->SetParameterName("A18Pos_y",false);
  a18yPosCmd->SetRange("A18Pos_y>=0. && A18Pos_y<15.0");
  a18yPosCmd->SetDefaultUnit( "cm" );
  a18yPosCmd->SetUnitCategory("Length");
  a18yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a19yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A19position",this);
  a19yPosCmd->SetGuidance("Set A19 position ");
  a19yPosCmd->SetParameterName("A19Pos_y",false);
  a19yPosCmd->SetRange("A19Pos_y>=0. && A19Pos_y<15.0");
  a19yPosCmd->SetDefaultUnit( "cm" );
  a19yPosCmd->SetUnitCategory("Length");
  a19yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a20yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A20position",this);
  a20yPosCmd->SetGuidance("Set A20 position ");
  a20yPosCmd->SetParameterName("A20Pos_y",false);
  a20yPosCmd->SetRange("A20Pos_y>=0. && A20Pos_y<15.0");
  a20yPosCmd->SetDefaultUnit( "cm" );
  a20yPosCmd->SetUnitCategory("Length");
  a20yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a21yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A21position",this);
  a21yPosCmd->SetGuidance("Set A21 position ");
  a21yPosCmd->SetParameterName("A21Pos_y",false);
  a21yPosCmd->SetRange("A21Pos_y>=0. && A21Pos_y<15.0");
  a21yPosCmd->SetDefaultUnit( "cm" );
  a21yPosCmd->SetUnitCategory("Length");
  a21yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a22yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A22position",this);
  a22yPosCmd->SetGuidance("Set A22 position ");
  a22yPosCmd->SetParameterName("A22Pos_y",false);
  a22yPosCmd->SetRange("A22Pos_y>=0. && A22Pos_y<15.0");
  a22yPosCmd->SetDefaultUnit( "cm" );
  a22yPosCmd->SetUnitCategory("Length");
  a22yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  a23yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A23position",this);
  a23yPosCmd->SetGuidance("Set A23 position ");
  a23yPosCmd->SetParameterName("A23Pos_y",false);
  a23yPosCmd->SetRange("A23Pos_y>=0. && A23Pos_y<15.0");
  a23yPosCmd->SetDefaultUnit( "cm" );
  a23yPosCmd->SetUnitCategory("Length");
  a23yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a24yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A24position",this);
  a24yPosCmd->SetGuidance("Set A24 position ");
  a24yPosCmd->SetParameterName("A24Pos_y",false);
  a24yPosCmd->SetRange("A24Pos_y>=0. && A24Pos_y<15.0");
  a24yPosCmd->SetDefaultUnit( "cm" );
  a24yPosCmd->SetUnitCategory("Length");
  a24yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a25yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A25position",this);
  a25yPosCmd->SetGuidance("Set A25 position ");
  a25yPosCmd->SetParameterName("A25Pos_y",false);
  a25yPosCmd->SetRange("A25Pos_y>=0. && A25Pos_y<15.0");
  a25yPosCmd->SetDefaultUnit( "cm" );
  a25yPosCmd->SetUnitCategory("Length");
  a25yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a26yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A26position",this);
  a26yPosCmd->SetGuidance("Set A26 position ");
  a26yPosCmd->SetParameterName("A26Pos_y",false);
  a26yPosCmd->SetRange("A26Pos_y>=0. && A26Pos_y<15.0");
  a26yPosCmd->SetDefaultUnit( "cm" );
  a26yPosCmd->SetUnitCategory("Length");
  a26yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a27yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A27position",this);
  a27yPosCmd->SetGuidance("Set A27 position ");
  a27yPosCmd->SetParameterName("A27Pos_y",false);
  a27yPosCmd->SetRange("A27Pos_y>=0. && A27Pos_y<15.0");
  a27yPosCmd->SetDefaultUnit( "cm" );
  a27yPosCmd->SetUnitCategory("Length");
  a27yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a28yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A28position",this);
  a28yPosCmd->SetGuidance("Set A28 position ");
  a28yPosCmd->SetParameterName("A28Pos_y",false);
  a28yPosCmd->SetRange("A28Pos_y>=0. && A28Pos_y<15.0");
  a28yPosCmd->SetDefaultUnit( "cm" );
  a28yPosCmd->SetUnitCategory("Length");
  a28yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a29yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A29position",this);
  a29yPosCmd->SetGuidance("Set A29 position ");
  a29yPosCmd->SetParameterName("A29Pos_y",false);
  a29yPosCmd->SetRange("A29Pos_y>=0. && A29Pos_y<15.0");
  a29yPosCmd->SetDefaultUnit( "cm" );
  a29yPosCmd->SetUnitCategory("Length");
  a29yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a30yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A30position",this);
  a30yPosCmd->SetGuidance("Set A30 position ");
  a30yPosCmd->SetParameterName("A30Pos_y",false);
  a30yPosCmd->SetRange("A30Pos_y>=0. && A30Pos_y<15.0");
  a30yPosCmd->SetDefaultUnit( "cm" );
  a30yPosCmd->SetUnitCategory("Length");
  a30yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  a31yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A31position",this);
  a31yPosCmd->SetGuidance("Set A31 position ");
  a31yPosCmd->SetParameterName("A31Pos_y",false);
  a31yPosCmd->SetRange("A31Pos_y>=0. && A31Pos_y<15.0");
  a31yPosCmd->SetDefaultUnit( "cm" );
  a31yPosCmd->SetUnitCategory("Length");
  a31yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a32yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A32position",this);
  a32yPosCmd->SetGuidance("Set A32 position ");
  a32yPosCmd->SetParameterName("A32Pos_y",false);
  a32yPosCmd->SetRange("A32Pos_y>=0. && A32Pos_y<15.0");
  a32yPosCmd->SetDefaultUnit( "cm" );
  a32yPosCmd->SetUnitCategory("Length");
  a32yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  a33yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A33position",this);
  a33yPosCmd->SetGuidance("Set A33 position ");
  a33yPosCmd->SetParameterName("A33Pos_y",false);
  a33yPosCmd->SetRange("A33Pos_y>=0. && A33Pos_y<15.0");
  a33yPosCmd->SetDefaultUnit( "cm" );
  a33yPosCmd->SetUnitCategory("Length");
  a33yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a34yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A34position",this);
  a34yPosCmd->SetGuidance("Set A34 position ");
  a34yPosCmd->SetParameterName("A34Pos_y",false);
  a34yPosCmd->SetRange("A34Pos_y>=0. && A34Pos_y<15.0");
  a34yPosCmd->SetDefaultUnit( "cm" );
  a34yPosCmd->SetUnitCategory("Length");
  a34yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a35yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A35position",this);
  a35yPosCmd->SetGuidance("Set A35 position ");
  a35yPosCmd->SetParameterName("A35Pos_y",false);
  a35yPosCmd->SetRange("A35Pos_y>=0. && A35Pos_y<15.0");
  a35yPosCmd->SetDefaultUnit( "cm" );
  a35yPosCmd->SetUnitCategory("Length");
  a35yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a36yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A36position",this);
  a36yPosCmd->SetGuidance("Set A36 position ");
  a36yPosCmd->SetParameterName("A36Pos_y",false);
  a36yPosCmd->SetRange("A36Pos_y>=0. && A36Pos_y<15.0");
  a36yPosCmd->SetDefaultUnit( "cm" );
  a36yPosCmd->SetUnitCategory("Length");
  a36yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a37yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A37position",this);
  a37yPosCmd->SetGuidance("Set A37 position ");
  a37yPosCmd->SetParameterName("A37Pos_y",false);
  a37yPosCmd->SetRange("A37Pos_y>=0. && A37Pos_y<15.0");
  a37yPosCmd->SetDefaultUnit( "cm" );
  a37yPosCmd->SetUnitCategory("Length");
  a37yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a38yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A38position",this);
  a38yPosCmd->SetGuidance("Set A38 position ");
  a38yPosCmd->SetParameterName("A38Pos_y",false);
  a38yPosCmd->SetRange("A38Pos_y>=0. && A38Pos_y<15.0");
  a38yPosCmd->SetDefaultUnit( "cm" );
  a38yPosCmd->SetUnitCategory("Length");
  a38yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  a39yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A39position",this);
  a39yPosCmd->SetGuidance("Set A39 position ");
  a39yPosCmd->SetParameterName("A39Pos_y",false);
  a39yPosCmd->SetRange("A39Pos_y>=0. && A39Pos_y<15.0");
  a39yPosCmd->SetDefaultUnit( "cm" );
  a39yPosCmd->SetUnitCategory("Length");
  a39yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  a40yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/A/A40position",this);
  a40yPosCmd->SetGuidance("Set A40 position ");
  a40yPosCmd->SetParameterName("A40Pos_y",false);
  a40yPosCmd->SetRange("A40Pos_y>=0. && A40Pos_y<15.0");
  a40yPosCmd->SetDefaultUnit( "cm" );
  a40yPosCmd->SetUnitCategory("Length");
  a40yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);





  b1yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B1position",this);
  b1yPosCmd->SetGuidance("Set B1 position ");
  b1yPosCmd->SetParameterName("B1Pos_y",false);
  b1yPosCmd->SetRange("B1Pos_y>=0. && B1Pos_y<15.0");
  b1yPosCmd->SetDefaultUnit( "cm" );
  b1yPosCmd->SetUnitCategory("Length");
  b1yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 

  b2yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B2position",this);
  b2yPosCmd->SetGuidance("Set B2 position ");
  b2yPosCmd->SetParameterName("B2Pos_y",false);
  b2yPosCmd->SetRange("B2Pos_y>=0. && B2Pos_y<15.0");
  b2yPosCmd->SetDefaultUnit( "cm" );
  b2yPosCmd->SetUnitCategory("Length");
  b2yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  b3yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B3position",this);
  b3yPosCmd->SetGuidance("Set B3 position ");
  b3yPosCmd->SetParameterName("B3Pos_y",false);
  b3yPosCmd->SetRange("B3Pos_y>=0. && B3Pos_y<15.0");
  b3yPosCmd->SetDefaultUnit( "cm" );
  b3yPosCmd->SetUnitCategory("Length");
  b3yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b4yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B4position",this);
  b4yPosCmd->SetGuidance("Set B4 position ");
  b4yPosCmd->SetParameterName("B4Pos_y",false);
  b4yPosCmd->SetRange("B4Pos_y>=0. && B4Pos_y<15.0");
  b4yPosCmd->SetDefaultUnit( "cm" );
  b4yPosCmd->SetUnitCategory("Length");
  b4yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b5yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B5position",this);
  b5yPosCmd->SetGuidance("Set B5 position ");
  b5yPosCmd->SetParameterName("B5Pos_y",false);
  b5yPosCmd->SetRange("B5Pos_y>=0. && B5Pos_y<15.0");
  b5yPosCmd->SetDefaultUnit( "cm" );
  b5yPosCmd->SetUnitCategory("Length");
  b5yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b6yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B6position",this);
  b6yPosCmd->SetGuidance("Set B6 position ");
  b6yPosCmd->SetParameterName("B6Pos_y",false);
  b6yPosCmd->SetRange("B6Pos_y>=0. && B6Pos_y<15.0");
  b6yPosCmd->SetDefaultUnit( "cm" );
  b6yPosCmd->SetUnitCategory("Length");
  b6yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b7yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B7position",this);
  b7yPosCmd->SetGuidance("Set B7 position ");
  b7yPosCmd->SetParameterName("B7Pos_y",false);
  b7yPosCmd->SetRange("B7Pos_y>=0. && B7Pos_y<15.0");
  b7yPosCmd->SetDefaultUnit( "cm" );
  b7yPosCmd->SetUnitCategory("Length");
  b7yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b8yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B8position",this);
  b8yPosCmd->SetGuidance("Set B8 position ");
  b8yPosCmd->SetParameterName("B8Pos_y",false);
  b8yPosCmd->SetRange("B8Pos_y>=0. && B8Pos_y<15.0");
  b8yPosCmd->SetDefaultUnit( "cm" );
  b8yPosCmd->SetUnitCategory("Length");
  b8yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b9yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B9position",this);
  b9yPosCmd->SetGuidance("Set B9 position ");
  b9yPosCmd->SetParameterName("B9Pos_y",false);
  b9yPosCmd->SetRange("B9Pos_y>=0. && B9Pos_y<15.0");
  b9yPosCmd->SetDefaultUnit( "cm" );
  b9yPosCmd->SetUnitCategory("Length");
  b9yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b10yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B10position",this);
  b10yPosCmd->SetGuidance("Set B10 position ");
  b10yPosCmd->SetParameterName("B10Pos_y",false);
  b10yPosCmd->SetRange("B10Pos_y>=0. && B10Pos_y<15.0");
  b10yPosCmd->SetDefaultUnit( "cm" );
  b10yPosCmd->SetUnitCategory("Length");
  b10yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b11yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B11position",this);
  b11yPosCmd->SetGuidance("Set B11 position ");
  b11yPosCmd->SetParameterName("B11Pos_y",false);
  b11yPosCmd->SetRange("B11Pos_y>=0. && B11Pos_y<15.0");
  b11yPosCmd->SetDefaultUnit( "cm" );
  b11yPosCmd->SetUnitCategory("Length");
  b11yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b12yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B12position",this);
  b12yPosCmd->SetGuidance("Set B12 position ");
  b12yPosCmd->SetParameterName("B12Pos_y",false);
  b12yPosCmd->SetRange("B12Pos_y>=0. && B12Pos_y<15.0");
  b12yPosCmd->SetDefaultUnit( "cm" );
  b12yPosCmd->SetUnitCategory("Length");
  b12yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b13yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B13position",this);
  b13yPosCmd->SetGuidance("Set B13 position ");
  b13yPosCmd->SetParameterName("B13Pos_y",false);
  b13yPosCmd->SetRange("B13Pos_y>=0. && B13Pos_y<15.0");
  b13yPosCmd->SetDefaultUnit( "cm" );
  b13yPosCmd->SetUnitCategory("Length");
  b13yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b14yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B14position",this);
  b14yPosCmd->SetGuidance("Set B14 position ");
  b14yPosCmd->SetParameterName("B14Pos_y",false);
  b14yPosCmd->SetRange("B14Pos_y>=0. && B14Pos_y<15.0");
  b14yPosCmd->SetDefaultUnit( "cm" );
  b14yPosCmd->SetUnitCategory("Length");
  b14yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b15yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B15position",this);
  b15yPosCmd->SetGuidance("Set B15 position ");
  b15yPosCmd->SetParameterName("B15Pos_y",false);
  b15yPosCmd->SetRange("B15Pos_y>=0. && B15Pos_y<15.0");
  b15yPosCmd->SetDefaultUnit( "cm" );
  b15yPosCmd->SetUnitCategory("Length");
  b15yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b16yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B16position",this);
  b16yPosCmd->SetGuidance("Set B16 position ");
  b16yPosCmd->SetParameterName("B16Pos_y",false);
  b16yPosCmd->SetRange("B16Pos_y>=0. && B16Pos_y<15.0");
  b16yPosCmd->SetDefaultUnit( "cm" );
  b16yPosCmd->SetUnitCategory("Length");
  b16yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b17yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B17position",this);
  b17yPosCmd->SetGuidance("Set B17 position ");
  b17yPosCmd->SetParameterName("B17Pos_y",false);
  b17yPosCmd->SetRange("B17Pos_y>=0. && B17Pos_y<15.0");
  b17yPosCmd->SetDefaultUnit( "cm" );
  b17yPosCmd->SetUnitCategory("Length");
  b17yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b18yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B18position",this);
  b18yPosCmd->SetGuidance("Set B18 position ");
  b18yPosCmd->SetParameterName("B18Pos_y",false);
  b18yPosCmd->SetRange("B18Pos_y>=0. && B18Pos_y<15.0");
  b18yPosCmd->SetDefaultUnit( "cm" );
  b18yPosCmd->SetUnitCategory("Length");
  b18yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b19yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B19position",this);
  b19yPosCmd->SetGuidance("Set B19 position ");
  b19yPosCmd->SetParameterName("B19Pos_y",false);
  b19yPosCmd->SetRange("B19Pos_y>=0. && B19Pos_y<15.0");
  b19yPosCmd->SetDefaultUnit( "cm" );
  b19yPosCmd->SetUnitCategory("Length");
  b19yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b20yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B20position",this);
  b20yPosCmd->SetGuidance("Set B20 position ");
  b20yPosCmd->SetParameterName("B20Pos_y",false);
  b20yPosCmd->SetRange("B20Pos_y>=0. && B20Pos_y<15.0");
  b20yPosCmd->SetDefaultUnit( "cm" );
  b20yPosCmd->SetUnitCategory("Length");
  b20yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b21yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B21position",this);
  b21yPosCmd->SetGuidance("Set B21 position ");
  b21yPosCmd->SetParameterName("B21Pos_y",false);
  b21yPosCmd->SetRange("B21Pos_y>=0. && B21Pos_y<15.0");
  b21yPosCmd->SetDefaultUnit( "cm" );
  b21yPosCmd->SetUnitCategory("Length");
  b21yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b22yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B22position",this);
  b22yPosCmd->SetGuidance("Set B22 position ");
  b22yPosCmd->SetParameterName("B22Pos_y",false);
  b22yPosCmd->SetRange("B22Pos_y>=0. && B22Pos_y<15.0");
  b22yPosCmd->SetDefaultUnit( "cm" );
  b22yPosCmd->SetUnitCategory("Length");
  b22yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  b23yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B23position",this);
  b23yPosCmd->SetGuidance("Set B23 position ");
  b23yPosCmd->SetParameterName("B23Pos_y",false);
  b23yPosCmd->SetRange("B23Pos_y>=0. && B23Pos_y<15.0");
  b23yPosCmd->SetDefaultUnit( "cm" );
  b23yPosCmd->SetUnitCategory("Length");
  b23yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  b24yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B24position",this);
  b24yPosCmd->SetGuidance("Set B24 position ");
  b24yPosCmd->SetParameterName("B24Pos_y",false);
  b24yPosCmd->SetRange("B24Pos_y>=0. && B24Pos_y<15.0");
  b24yPosCmd->SetDefaultUnit( "cm" );
  b24yPosCmd->SetUnitCategory("Length");
  b24yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b25yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B25position",this);
  b25yPosCmd->SetGuidance("Set B25 position ");
  b25yPosCmd->SetParameterName("B25Pos_y",false);
  b25yPosCmd->SetRange("B25Pos_y>=0. && B25Pos_y<15.0");
  b25yPosCmd->SetDefaultUnit( "cm" );
  b25yPosCmd->SetUnitCategory("Length");
  b25yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 

  b26yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B26position",this);
  b26yPosCmd->SetGuidance("Set B26 position ");
  b26yPosCmd->SetParameterName("B26Pos_y",false);
  b26yPosCmd->SetRange("B26Pos_y>=0. && B26Pos_y<15.0");
  b26yPosCmd->SetDefaultUnit( "cm" );
  b26yPosCmd->SetUnitCategory("Length");
  b26yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b27yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B27position",this);
  b27yPosCmd->SetGuidance("Set B27 position ");
  b27yPosCmd->SetParameterName("B27Pos_y",false);
  b27yPosCmd->SetRange("B27Pos_y>=0. && B27Pos_y<15.0");
  b27yPosCmd->SetDefaultUnit( "cm" );
  b27yPosCmd->SetUnitCategory("Length");
  b27yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b28yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B28position",this);
  b28yPosCmd->SetGuidance("Set B28 position ");
  b28yPosCmd->SetParameterName("B28Pos_y",false);
  b28yPosCmd->SetRange("B28Pos_y>=0. && B28Pos_y<15.0");
  b28yPosCmd->SetDefaultUnit( "cm" );
  b28yPosCmd->SetUnitCategory("Length");
  b28yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b29yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B29position",this);
  b29yPosCmd->SetGuidance("Set B29 position ");
  b29yPosCmd->SetParameterName("B29Pos_y",false);
  b29yPosCmd->SetRange("B29Pos_y>=0. && B29Pos_y<15.0");
  b29yPosCmd->SetDefaultUnit( "cm" );
  b29yPosCmd->SetUnitCategory("Length");
  b29yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b30yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B30position",this);
  b30yPosCmd->SetGuidance("Set B30 position ");
  b30yPosCmd->SetParameterName("B30Pos_y",false);
  b30yPosCmd->SetRange("B30Pos_y>=0. && B30Pos_y<15.0");
  b30yPosCmd->SetDefaultUnit( "cm" );
  b30yPosCmd->SetUnitCategory("Length");
  b30yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  b31yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B31position",this);
  b31yPosCmd->SetGuidance("Set B31 position ");
  b31yPosCmd->SetParameterName("B31Pos_y",false);
  b31yPosCmd->SetRange("B31Pos_y>=0. && B31Pos_y<15.0");
  b31yPosCmd->SetDefaultUnit( "cm" );
  b31yPosCmd->SetUnitCategory("Length");
  b31yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b32yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B32position",this);
  b32yPosCmd->SetGuidance("Set B32 position ");
  b32yPosCmd->SetParameterName("B32Pos_y",false);
  b32yPosCmd->SetRange("B32Pos_y>=0. && B32Pos_y<15.0");
  b32yPosCmd->SetDefaultUnit( "cm" );
  b32yPosCmd->SetUnitCategory("Length");
  b32yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b33yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B33position",this);
  b33yPosCmd->SetGuidance("Set B33 position ");
  b33yPosCmd->SetParameterName("B33Pos_y",false);
  b33yPosCmd->SetRange("B33Pos_y>=0. && B33Pos_y<15.0");
  b33yPosCmd->SetDefaultUnit( "cm" );
  b33yPosCmd->SetUnitCategory("Length");
  b33yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b34yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B34position",this);
  b34yPosCmd->SetGuidance("Set B34 position ");
  b34yPosCmd->SetParameterName("B34Pos_y",false);
  b34yPosCmd->SetRange("B34Pos_y>=0. && B34Pos_y<15.0");
  b34yPosCmd->SetDefaultUnit( "cm" );
  b34yPosCmd->SetUnitCategory("Length");
  b34yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b35yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B35position",this);
  b35yPosCmd->SetGuidance("Set B35 position ");
  b35yPosCmd->SetParameterName("B35Pos_y",false);
  b35yPosCmd->SetRange("B35Pos_y>=0. && B35Pos_y<15.0");
  b35yPosCmd->SetDefaultUnit( "cm" );
  b35yPosCmd->SetUnitCategory("Length");
  b35yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b36yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B36position",this);
  b36yPosCmd->SetGuidance("Set B36 position ");
  b36yPosCmd->SetParameterName("B36Pos_y",false);
  b36yPosCmd->SetRange("B36Pos_y>=0. && B36Pos_y<15.0");
  b36yPosCmd->SetDefaultUnit( "cm" );
  b36yPosCmd->SetUnitCategory("Length");
  b36yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b37yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B37position",this);
  b37yPosCmd->SetGuidance("Set B37 position ");
  b37yPosCmd->SetParameterName("B37Pos_y",false);
  b37yPosCmd->SetRange("B37Pos_y>=0. && B37Pos_y<15.0");
  b37yPosCmd->SetDefaultUnit( "cm" );
  b37yPosCmd->SetUnitCategory("Length");
  b37yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b38yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B38position",this);
  b38yPosCmd->SetGuidance("Set B38 position ");
  b38yPosCmd->SetParameterName("B38Pos_y",false);
  b38yPosCmd->SetRange("B38Pos_y>=0. && B38Pos_y<15.0");
  b38yPosCmd->SetDefaultUnit( "cm" );
  b38yPosCmd->SetUnitCategory("Length");
  b38yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  b39yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B39position",this);
  b39yPosCmd->SetGuidance("Set B39 position ");
  b39yPosCmd->SetParameterName("B39Pos_y",false);
  b39yPosCmd->SetRange("B39Pos_y>=0. && B39Pos_y<15.0");
  b39yPosCmd->SetDefaultUnit( "cm" );
  b39yPosCmd->SetUnitCategory("Length");
  b39yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  b40yPosCmd = new G4UIcmdWithADoubleAndUnit("/MLC/B/B40position",this);
  b40yPosCmd->SetGuidance("Set B40 position ");
  b40yPosCmd->SetParameterName("B40Pos_y",false);
  b40yPosCmd->SetRange("B40Pos_y>=0. && B40Pos_y<15.0");
  b40yPosCmd->SetDefaultUnit( "cm" );
  b40yPosCmd->SetUnitCategory("Length");
  b40yPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//****************************************************************************
MedLinacMLCMessenger::~MedLinacMLCMessenger()
{
  delete a1yPosCmd;
  delete a2yPosCmd;
  delete a3yPosCmd;
  delete a4yPosCmd;
  delete a5yPosCmd;
  delete a6yPosCmd;
  delete a7yPosCmd;
  delete a8yPosCmd;
  delete a9yPosCmd;
  delete a10yPosCmd;
  delete a11yPosCmd;
  delete a12yPosCmd;
  delete a13yPosCmd;
  delete a14yPosCmd;
  delete a15yPosCmd;
  delete a16yPosCmd;
  delete a17yPosCmd;
  delete a18yPosCmd;
  delete a19yPosCmd;
  delete a20yPosCmd;
  delete a21yPosCmd;
  delete a22yPosCmd;
  delete a23yPosCmd;
  delete a24yPosCmd;
  delete a25yPosCmd;
  delete a26yPosCmd;
  delete a27yPosCmd;
  delete a28yPosCmd;
  delete a29yPosCmd;
  delete a30yPosCmd;
  delete a31yPosCmd;
  delete a32yPosCmd;
  delete a33yPosCmd;
  delete a34yPosCmd;
  delete a35yPosCmd;
  delete a36yPosCmd;
  delete a37yPosCmd;
  delete a38yPosCmd;
  delete a39yPosCmd;
  delete a40yPosCmd;

  delete b1yPosCmd;
  delete b2yPosCmd;
  delete b3yPosCmd;
  delete b4yPosCmd;
  delete b5yPosCmd;
  delete b6yPosCmd;
  delete b7yPosCmd;
  delete b8yPosCmd;
  delete b9yPosCmd;
  delete b10yPosCmd;
  delete b11yPosCmd;
  delete b12yPosCmd;
  delete b13yPosCmd;
  delete b14yPosCmd;
  delete b15yPosCmd;
  delete b16yPosCmd;
  delete b17yPosCmd;
  delete b18yPosCmd;
  delete b19yPosCmd;
  delete b20yPosCmd;
  delete b21yPosCmd;
  delete b22yPosCmd;
  delete b23yPosCmd;
  delete b24yPosCmd;
  delete b25yPosCmd;
  delete b26yPosCmd;
  delete b27yPosCmd;
  delete b28yPosCmd;
  delete b29yPosCmd;
  delete b30yPosCmd;
  delete b31yPosCmd;
  delete b32yPosCmd;
  delete b33yPosCmd;
  delete b34yPosCmd;
  delete b35yPosCmd;
  delete b36yPosCmd;
  delete b37yPosCmd;
  delete b38yPosCmd;
  delete b39yPosCmd;
  delete b40yPosCmd;

  delete MedLinacMLCDir;
  
  delete aDir;
}

//****************************************************************************

void MedLinacMLCMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == a1yPosCmd )
   { MedLinacMLCDeco->SetA1Pos_y(a1yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a2yPosCmd )
   { MedLinacMLCDeco->SetA2Pos_y(a2yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a3yPosCmd )
   { MedLinacMLCDeco->SetA3Pos_y(a3yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a4yPosCmd )
   { MedLinacMLCDeco->SetA4Pos_y(a4yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a5yPosCmd )
   { MedLinacMLCDeco->SetA5Pos_y(a5yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a6yPosCmd )
   { MedLinacMLCDeco->SetA6Pos_y(a6yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a7yPosCmd )
   { MedLinacMLCDeco->SetA7Pos_y(a7yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a8yPosCmd )
   { MedLinacMLCDeco->SetA8Pos_y(a8yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a9yPosCmd )
   { MedLinacMLCDeco->SetA9Pos_y(a9yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a10yPosCmd )
   { MedLinacMLCDeco->SetA10Pos_y(a10yPosCmd->GetNewDoubleValue(newValue));}

  if( command == a11yPosCmd )
   { MedLinacMLCDeco->SetA11Pos_y(a11yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a12yPosCmd )
   { MedLinacMLCDeco->SetA12Pos_y(a12yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a13yPosCmd )
   { MedLinacMLCDeco->SetA13Pos_y(a13yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a14yPosCmd )
   { MedLinacMLCDeco->SetA14Pos_y(a14yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a15yPosCmd )
   { MedLinacMLCDeco->SetA15Pos_y(a15yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a16yPosCmd )
   { MedLinacMLCDeco->SetA16Pos_y(a16yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a17yPosCmd )
   { MedLinacMLCDeco->SetA17Pos_y(a17yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a18yPosCmd )
   { MedLinacMLCDeco->SetA18Pos_y(a18yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a19yPosCmd )
   { MedLinacMLCDeco->SetA19Pos_y(a19yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a20yPosCmd )
   { MedLinacMLCDeco->SetA20Pos_y(a20yPosCmd->GetNewDoubleValue(newValue));}

  if( command == a21yPosCmd )
   { MedLinacMLCDeco->SetA21Pos_y(a21yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a22yPosCmd )
   { MedLinacMLCDeco->SetA22Pos_y(a22yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a23yPosCmd )
   { MedLinacMLCDeco->SetA23Pos_y(a23yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a24yPosCmd )
   { MedLinacMLCDeco->SetA24Pos_y(a24yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a25yPosCmd )
   { MedLinacMLCDeco->SetA25Pos_y(a25yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a26yPosCmd )
   { MedLinacMLCDeco->SetA26Pos_y(a26yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a27yPosCmd )
   { MedLinacMLCDeco->SetA27Pos_y(a27yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a28yPosCmd )
   { MedLinacMLCDeco->SetA28Pos_y(a28yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a29yPosCmd )
   { MedLinacMLCDeco->SetA29Pos_y(a29yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a30yPosCmd )
   { MedLinacMLCDeco->SetA30Pos_y(a30yPosCmd->GetNewDoubleValue(newValue));}

  if( command == a31yPosCmd )
   { MedLinacMLCDeco->SetA31Pos_y(a31yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a32yPosCmd )
   { MedLinacMLCDeco->SetA32Pos_y(a32yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a33yPosCmd )
   { MedLinacMLCDeco->SetA33Pos_y(a33yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a34yPosCmd )
   { MedLinacMLCDeco->SetA34Pos_y(a34yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a35yPosCmd )
   { MedLinacMLCDeco->SetA35Pos_y(a35yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a36yPosCmd )
   { MedLinacMLCDeco->SetA36Pos_y(a36yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a37yPosCmd )
   { MedLinacMLCDeco->SetA37Pos_y(a37yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a38yPosCmd )
   { MedLinacMLCDeco->SetA38Pos_y(a38yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a39yPosCmd )
   { MedLinacMLCDeco->SetA39Pos_y(a39yPosCmd->GetNewDoubleValue(newValue));}
  if( command == a40yPosCmd )
   { MedLinacMLCDeco->SetA40Pos_y(a40yPosCmd->GetNewDoubleValue(newValue));}







  if( command == b1yPosCmd )
   { MedLinacMLCDeco->SetB1Pos_y(b1yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b2yPosCmd )
   { MedLinacMLCDeco->SetB2Pos_y(b2yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b3yPosCmd )
   { MedLinacMLCDeco->SetB3Pos_y(b3yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b4yPosCmd )
   { MedLinacMLCDeco->SetB4Pos_y(b4yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b5yPosCmd )
   { MedLinacMLCDeco->SetB5Pos_y(b5yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b6yPosCmd )
   { MedLinacMLCDeco->SetB6Pos_y(b6yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b7yPosCmd )
   { MedLinacMLCDeco->SetB7Pos_y(b7yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b8yPosCmd )
   { MedLinacMLCDeco->SetB8Pos_y(b8yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b9yPosCmd )
   { MedLinacMLCDeco->SetB9Pos_y(b9yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b10yPosCmd )
   { MedLinacMLCDeco->SetB10Pos_y(b10yPosCmd->GetNewDoubleValue(newValue));}

  if( command == b11yPosCmd )
   { MedLinacMLCDeco->SetB11Pos_y(b11yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b12yPosCmd )
   { MedLinacMLCDeco->SetB12Pos_y(b12yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b13yPosCmd )
   { MedLinacMLCDeco->SetB13Pos_y(b13yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b14yPosCmd )
   { MedLinacMLCDeco->SetB14Pos_y(b14yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b15yPosCmd )
   { MedLinacMLCDeco->SetB15Pos_y(b15yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b16yPosCmd )
   { MedLinacMLCDeco->SetB16Pos_y(b16yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b17yPosCmd )
   { MedLinacMLCDeco->SetB17Pos_y(b17yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b18yPosCmd )
   { MedLinacMLCDeco->SetB18Pos_y(b18yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b19yPosCmd )
   { MedLinacMLCDeco->SetB19Pos_y(b19yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b20yPosCmd )
   { MedLinacMLCDeco->SetB20Pos_y(b20yPosCmd->GetNewDoubleValue(newValue));}

  if( command == b21yPosCmd )
   { MedLinacMLCDeco->SetB21Pos_y(b21yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b22yPosCmd )
   { MedLinacMLCDeco->SetB22Pos_y(b22yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b23yPosCmd )
   { MedLinacMLCDeco->SetB23Pos_y(b23yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b24yPosCmd )
   { MedLinacMLCDeco->SetB24Pos_y(b24yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b25yPosCmd )
   { MedLinacMLCDeco->SetB25Pos_y(b25yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b26yPosCmd )
   { MedLinacMLCDeco->SetB26Pos_y(b26yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b27yPosCmd )
   { MedLinacMLCDeco->SetB27Pos_y(b27yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b28yPosCmd )
   { MedLinacMLCDeco->SetB28Pos_y(b28yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b29yPosCmd )
   { MedLinacMLCDeco->SetB29Pos_y(b29yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b30yPosCmd )
   { MedLinacMLCDeco->SetB30Pos_y(b30yPosCmd->GetNewDoubleValue(newValue));}

  if( command == b31yPosCmd )
   { MedLinacMLCDeco->SetB31Pos_y(b31yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b32yPosCmd )
   { MedLinacMLCDeco->SetB32Pos_y(b32yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b33yPosCmd )
   { MedLinacMLCDeco->SetB33Pos_y(b33yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b34yPosCmd )
   { MedLinacMLCDeco->SetB34Pos_y(b34yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b35yPosCmd )
   { MedLinacMLCDeco->SetB35Pos_y(b35yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b36yPosCmd )
   { MedLinacMLCDeco->SetB36Pos_y(b36yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b37yPosCmd )
   { MedLinacMLCDeco->SetB37Pos_y(b37yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b38yPosCmd )
   { MedLinacMLCDeco->SetB38Pos_y(b38yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b39yPosCmd )
   { MedLinacMLCDeco->SetB39Pos_y(b39yPosCmd->GetNewDoubleValue(newValue));}
  if( command == b40yPosCmd )
   { MedLinacMLCDeco->SetB40Pos_y(b40yPosCmd->GetNewDoubleValue(newValue));}
}

//****************************************************************************
