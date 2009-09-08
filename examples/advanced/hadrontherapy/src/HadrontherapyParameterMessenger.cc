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
// $Id: HadrontherapyInteractionParametersMessenger.cc;
//

#include "HadrontherapyParameterMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"

#include "G4UIdirectory.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
G4double E=1.;

HadrontherapyParameterMessenger::HadrontherapyParameterMessenger(HadrontherapyInteractionParameters* param)
:pParam(param)
{
	paramDir = new G4UIdirectory("/parameter/");
    paramDir -> SetGuidance("commands to generate stopping power and range");
    
	dedxCmd = new G4UIcmdWithAString("/parameter/getstopping",this);  
    dedxCmd->SetGuidance("Get mass stopping power\n parameters: material minimum and maximum kinetic energy, n points, particle, filename ");
    dedxCmd->SetParameterName("inputData",false);
    dedxCmd->AvailableForStates(G4State_Idle);  

	rangeCmd = new G4UIcmdWithAString("/parameter/getrange",this);  
    rangeCmd->SetGuidance("Get CSDA Range\n parameters: minimum and maximum kinetic energy, n points, particle, material, filename ");
    rangeCmd->SetParameterName("inputData",false);
    rangeCmd->AvailableForStates(G4State_Idle);  

	materCmd = new G4UIcmdWithoutParameter("/parameter/materials",this);  
    materCmd->SetGuidance("Print Nist materials list");
    //materCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle,G4State_GeomClosed, G4State_EventProc);  
    materCmd->AvailableForStates(G4State_Idle);  
}
HadrontherapyParameterMessenger::~HadrontherapyParameterMessenger()
{
	delete paramDir;
	delete dedxCmd;
	delete rangeCmd;
	delete materCmd;

}

void HadrontherapyParameterMessenger::SetNewValue(G4UIcommand* command, G4String vararg)
{
	if (command == dedxCmd)
	{
		pParam -> GetStoppingTable(vararg);
	}
	else if (command == rangeCmd)
	{
		pParam -> GetCSDARangeTable(vararg);
	}
	else if (command == materCmd)
	{
		pParam -> ListOfNistMaterials();
	}

}

