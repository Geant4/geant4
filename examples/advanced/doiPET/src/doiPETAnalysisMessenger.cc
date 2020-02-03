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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan


#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "doiPETAnalysis.hh"
#include "doiPETAnalysisMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"


////////// Constructor //////////////////////////////////////////////
doiPETAnalysisMessenger::doiPETAnalysisMessenger(doiPETAnalysis* analy)
	:ptrAnalysis (analy)
{


	Dir = new G4UIdirectory("/change/");
	Dir->SetGuidance("Change parameters");

	//Change activity in the phantom
	changeActivityCmd = new G4UIcmdWithADoubleAndUnit("/change/Activity",this);
	changeActivityCmd->SetGuidance("Set activity in the phantom");
	changeActivityCmd->SetParameterName("activity",false);
	changeActivityCmd->SetUnitCategory("Activity");
	changeActivityCmd->SetRange("activity>0");
	//changeActivityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


	//Change Halflife of the isotope.
	changeHalfLifeCmd = new G4UIcmdWithADoubleAndUnit("/change/HalfLife",this);
	changeHalfLifeCmd->SetGuidance("Set Halflife of the isotope");
	changeHalfLifeCmd->SetParameterName("Time",false);
	changeHalfLifeCmd->SetUnitCategory("Time");

	//change the type of output: singles or coincidence out put
	//ChooseOutputTypeCmd = new G4UIcmdWithAString("/change/OutputType",this);
	//ChooseOutputTypeCmd->SetGuidance("Select the type of application");
	//ChooseOutputTypeCmd->SetParameterName("choice",false);
}

////////// Destructor //////////////////////////////////////////////
doiPETAnalysisMessenger::~doiPETAnalysisMessenger()
{	
	delete Dir;
	delete changeActivityCmd;
	delete changeHalfLifeCmd;
	//delete ChooseOutputTypeCmd;
}

////////// SetNewValue /////////////////////////////////////////////
void doiPETAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	/*if(command == ChooseOutputTypeCmd){
	ptrAnalysis->TypeOfOutput(newValue);
	}*/
	if (command == changeActivityCmd)
	{
		ptrAnalysis->SetActivity(changeActivityCmd->GetNewDoubleValue(newValue));
	}
	else if (command == changeHalfLifeCmd)
	{
		ptrAnalysis->SetIsotopeHalfLife(changeHalfLifeCmd->GetNewDoubleValue(newValue));
	}	
}


