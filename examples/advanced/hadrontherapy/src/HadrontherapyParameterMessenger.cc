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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyParameterMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "G4SystemOfUnits.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

HadrontherapyParameterMessenger::HadrontherapyParameterMessenger(HadrontherapyInteractionParameters* param)
:pParam(param)
{
    paramDir = new G4UIdirectory("/parameter/");
    paramDir -> SetGuidance("Commands to generate stopping power and range");
    
    dedxCmd = new G4UIcmdWithAString("/parameter/getstopping",this);  
    dedxCmd->SetGuidance("Get mass stopping powers"
			"\n[usage]: /parameter/getstopping Material [Emin] [Emax] [N] [Particle] [File]" 
			"\n         Material:(string) Material name, like G4_H, G4_WATER,..., look at /parameter/nist"
			"\n         Emin Emax:(double) minimum and maximum kinetic energy (MeV)"
			"\n         N:(double) [number of points]"
			"\n         Particle:(string) Particle name, look at /particle/list"
			"\n         File:(string) Name for the output file."
			"\nDefault values for parameters inside [] are respectively:"
			"\n \"1 MeV\", \"Emin\", \"1\", \"proton\", \"stdout\"");
    dedxCmd->SetParameterName("inputData",false);
    dedxCmd->AvailableForStates(G4State_Idle);  

    listCmd = new G4UIcmdWithAString("/parameter/nist",this);  
    listCmd -> SetGuidance("Print NIST elements/materials.\nParameters:"
			    "\n\t all: will print elements and compounds"
			    "\n\t simple: will print elements only"
			    "\n\t compound: will print compounds only"
			    "\n\t hep: will print hep compounds"
			    "\n\t list: will print a simple full list of all elements and compounds");
    listCmd -> SetParameterName("String",true);
    listCmd -> SetDefaultValue("list");
    listCmd -> SetCandidates("all simple compound hep list");
    listCmd ->AvailableForStates(G4State_Idle);  
    //Available G4 States (G4State_PreInit, G4State_Init, G4State_Idle,G4State_GeomClosed, G4State_EventProc);  
}
HadrontherapyParameterMessenger::~HadrontherapyParameterMessenger()
{
	delete paramDir;
	delete dedxCmd;
	delete listCmd;
}

void HadrontherapyParameterMessenger::SetNewValue(G4UIcommand* command, G4String vararg)
{
	if (command == dedxCmd)
	{
		pParam -> GetStoppingTable(vararg);
	}
	else if (command == listCmd)
	{
		pParam -> ListOfNistMaterials(vararg);
	}
}

