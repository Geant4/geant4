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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "DetectorMessenger.hh"

DetectorMessenger::~DetectorMessenger()
{
	// UI commands and directory have to be deleted
	delete fieldSide;
	delete sourceToSkinDistance;
	delete phantomSide;
	delete voxelSide;
	delete voxelDepth;
}

DetectorMessenger::DetectorMessenger(DetectorConstruction* myDetector)
: detPointer(myDetector)
{
	fieldSide = new G4UIcmdWithADoubleAndUnit("/DetectorConstruction/Acc/fieldSide", this);
	fieldSide -> SetDefaultUnit("mm");
	fieldSide -> SetDefaultValue(100.);
	fieldSide -> SetGuidance("Set the field side");
	detPointer -> SetJaws(100.*mm);

	sourceToSkinDistance = new G4UIcmdWithADoubleAndUnit("/DetectorConstruction/Acc/sourceToSkinDistance", this);
	sourceToSkinDistance -> SetDefaultUnit("mm");
	sourceToSkinDistance -> SetDefaultValue(900.);
	sourceToSkinDistance -> SetGuidance("Set the distance between the source and the patient's skin");
	detPointer -> SetTargetPosition(900.*mm);

	phantomSide = new G4UIcmdWithADoubleAndUnit("/DetectorConstruction/Phantom/phantomSide", this);
	phantomSide -> SetDefaultUnit("mm");
	phantomSide -> SetDefaultValue(255.*mm);
	phantomSide -> SetGuidance("Set the side of the phantom");
	detPointer -> SetPhantomSide(255.*mm);

	voxelSide = new G4UIcmdWithADoubleAndUnit("/DetectorConstruction/Phantom/voxelSide", this);
	voxelSide -> SetDefaultUnit("mm");
	voxelSide -> SetDefaultValue(5.*mm);
	voxelSide -> SetGuidance("Set the phantom voxel lateral dimension");
	detPointer -> SetVoxelSide(5.*mm);

	voxelDepth = new G4UIcmdWithADoubleAndUnit("/DetectorConstruction/Phantom/voxelDepth", this);
	voxelDepth -> SetDefaultUnit("mm");
	voxelDepth -> SetDefaultValue(5.*mm);
	voxelDepth -> SetGuidance("Set the phantom voxel depth");
	detPointer -> SetVoxelDepth(5.*mm);
}

void DetectorMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd == fieldSide)
	{
		fieldSide -> GetNewUnitValue(newValue);
		detPointer -> SetJaws(fieldSide -> GetNewDoubleValue(newValue));
		detPointer -> UpdateGeometry("fieldSide", fieldSide -> GetNewDoubleRawValue(newValue));
	}
	else if (cmd == sourceToSkinDistance)
	{
		sourceToSkinDistance -> GetNewUnitValue(newValue);
		detPointer -> SetTargetPosition(sourceToSkinDistance -> GetNewDoubleValue(newValue));
		detPointer -> UpdateGeometry("sourceToSkinDistance", sourceToSkinDistance -> GetNewDoubleRawValue(newValue));
	}
	else if (cmd == voxelDepth)
	{
		voxelDepth -> GetNewUnitValue(newValue);
		detPointer -> SetVoxelDepth(voxelDepth ->GetNewDoubleValue(newValue));
		detPointer -> UpdateGeometry("voxelDepth", voxelDepth -> GetNewDoubleRawValue(newValue));
	}
	else if (cmd == voxelSide)
	{
		voxelSide -> GetNewUnitValue(newValue);
		detPointer -> SetVoxelSide(voxelSide ->GetNewDoubleValue(newValue));
		detPointer -> UpdateGeometry("voxelSide", voxelSide -> GetNewDoubleRawValue(newValue));

	}
	else if (cmd == phantomSide)
	{
		phantomSide -> GetNewUnitValue(newValue);
		detPointer -> SetPhantomSide(phantomSide ->GetNewDoubleValue(newValue));
		detPointer -> UpdateGeometry("phantomSide", phantomSide -> GetNewDoubleRawValue(newValue));
	}
	else
		G4cerr << "DetectorMessenger::SetNewValue: command not found" << G4endl;
}




