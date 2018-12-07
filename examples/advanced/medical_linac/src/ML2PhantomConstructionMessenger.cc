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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2PhantomConstructionMessenger.hh"
#include "ML2PhantomConstruction.hh"

#include "G4ios.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


CML2PhantomConstructionMessenger::CML2PhantomConstructionMessenger(CML2PhantomConstruction *phantomConstructor) : pPhantomConstructor (phantomConstructor)
{
	PhantomName = new G4UIcmdWithAString("/phantom/PhantomName",this);
	PhantomName -> SetDefaultValue("fullWater");
	PhantomName -> SetGuidance("phantom name to select among those implemented (fullWater, boxInBox)");
	pPhantomConstructor -> setPhantomName("fullWater");

	phantomCentre = new G4UIcmdWith3VectorAndUnit("/phantom/centre", this);
	phantomCentre -> SetDefaultUnit("mm");
	phantomCentre -> SetGuidance("phantom centre coordinates in the world [mm]");
	phantomCentre -> SetDefaultValue(G4ThreeVector(0.,0.,0.));
}

CML2PhantomConstructionMessenger::~CML2PhantomConstructionMessenger(void)
{
	delete PhantomName;
	delete phantomCentre;
}
void CML2PhantomConstructionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd == PhantomName)
	{
		pPhantomConstructor->setPhantomName(newValue);
	}

	if (cmd == PhantomFileName )
	{
		pPhantomConstructor -> setPhantomFileName (newValue);
	}

	if (cmd == phantomCentre )
	{
	        pPhantomConstructor -> addNewCentre(phantomCentre->GetNew3VectorRawValue(newValue));
	}
}

