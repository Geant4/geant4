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
#include "AnalysisMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

AnalysisMessenger::AnalysisMessenger()
{
	analysisDir = new G4UIdirectory("/analysis/");
	analysisDir -> SetGuidance("Customize data analysis output");
	
	useRootForOutput = new G4UIcmdWithABool("/analysis/useRoot", this);
	useRootForOutput -> SetGuidance("If true, a ROOT file is used for the output. If false, a plaintext csv file is used");
	useRootForOutput -> SetParameterName("Root", true);
	useRootForOutput -> AvailableForStates(G4State_PreInit);
	
	//default
	rootOutput = true;
}

AnalysisMessenger::~AnalysisMessenger()
{
	delete useRootForOutput;
	
	delete analysisDir;
}

void AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{

	if( command == useRootForOutput )
	{
		rootOutput = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( rootOutput == true ) G4cout << "Outputting to a ROOT file" << G4endl;
		else G4cout << "Outputting to a csv file" << G4endl;
	}
}
