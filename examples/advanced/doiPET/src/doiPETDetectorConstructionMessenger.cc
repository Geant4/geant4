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

#include "doiPETDetectorConstruction.hh"
#include "doiPETDetectorConstructionMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"


////////// Constructor //////////////////////////////////////////////
doiPETDetectorConstructionMessenger::doiPETDetectorConstructionMessenger(doiPETDetectorConstruction* det)
	:Detector (det)
{
	Dir = new G4UIdirectory("/changePhantom/");
	Dir->SetGuidance("Choose the phantom");

	ChoosePhantomCmd = new G4UIcmdWithAString("/changePhantom/setPhantom",this);
	ChoosePhantomCmd->SetGuidance("Select Phantom type");
	ChoosePhantomCmd->SetParameterName("choice",false);
	ChoosePhantomCmd->AvailableForStates(G4State_PreInit);

	// Change Phantom position
	changeThePhantomPositionCmd = new G4UIcmdWith3VectorAndUnit("/changePhantom/position", this);
	changeThePhantomPositionCmd -> SetGuidance("Insert X Y and Z dimensions for the position of the center of the Phantom" " respect to that of the \"World\""); 
	changeThePhantomPositionCmd -> SetParameterName("PositionAlongX", "PositionAlongY", "PositionAlongZ", false);
	changeThePhantomPositionCmd -> SetDefaultUnit("mm");
	changeThePhantomPositionCmd -> SetUnitCandidates("um mm cm m"); 
	changeThePhantomPositionCmd -> AvailableForStates(G4State_PreInit);

	//Change Phantom Length and Radius
	changePhantomRadiusCmd = new G4UIcmdWithADoubleAndUnit("/changePhantom/Radius",this);
	changePhantomRadiusCmd->SetGuidance("Set radius of the phantom");
	changePhantomRadiusCmd->SetParameterName("radius",false);
	changePhantomRadiusCmd->SetUnitCategory("Length");
	changePhantomRadiusCmd->SetRange("radius>0");
	changePhantomRadiusCmd->AvailableForStates(G4State_PreInit);

	//Change the length of the phantom 
	changePhantomLengthCmd = new G4UIcmdWithADoubleAndUnit("/changePhantom/Length",this);
	changePhantomLengthCmd->SetGuidance("Set length of the phantom");
	changePhantomLengthCmd->SetParameterName("length",false);
	changePhantomLengthCmd->SetUnitCategory("Length");
	changePhantomLengthCmd->SetRange("length>0");
	changePhantomLengthCmd->AvailableForStates(G4State_PreInit);

}

////////// Destructor //////////////////////////////////////////////
doiPETDetectorConstructionMessenger::~doiPETDetectorConstructionMessenger()
{	
	delete Dir;
	delete ChoosePhantomCmd;
		
	delete changeThePhantomPositionCmd;
	delete changePhantomRadiusCmd;
	delete changePhantomLengthCmd;


}

////////// SetNewValue /////////////////////////////////////////////
void doiPETDetectorConstructionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	if ( command == ChoosePhantomCmd ) Detector->ChangePhantom(newValue);
	else if (command == changeThePhantomPositionCmd )
	{
		G4ThreeVector size = changeThePhantomPositionCmd -> GetNew3VectorValue(newValue);
		Detector -> SetPhantomPosition(size);
	}
	else if (command == changePhantomRadiusCmd)
	{
		Detector->SetPhantomRadius(changePhantomRadiusCmd->GetNewDoubleValue(newValue));
	}
	else if (command == changePhantomLengthCmd)
	{
		Detector->SetPhantomLength(changePhantomLengthCmd->GetNewDoubleValue(newValue));
	}
}


