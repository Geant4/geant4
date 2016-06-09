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
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2PhantomConstruction.hh"
#include "ML2PhantomConstructionMessenger.hh"
#include "G4UImanager.hh"

CML2PhantomConstruction::CML2PhantomConstruction(void)
{
	this->phantomContstructionMessenger=new CML2PhantomConstructionMessenger(this);
}

CML2PhantomConstruction::~CML2PhantomConstruction(void)
{
	if (this->phantomName=="fullWater")
	{
		delete Ph_fullWater;
	}
	else if (this->phantomName=="BoxInBox")
	{
		delete Ph_BoxInBox;
	}
}

CML2PhantomConstruction* CML2PhantomConstruction::instance = 0;

CML2PhantomConstruction* CML2PhantomConstruction::GetInstance(void)
{
  if (instance == 0)
    {
      instance = new CML2PhantomConstruction();
     
    }
  return instance;
}
void CML2PhantomConstruction::design(void)
{
// switch between two different phantoms according to the main macro selection
	if (this->phantomName=="fullWater")
	{
		this->Ph_fullWater=new CML2Ph_FullWater();
	}
	else if (this->phantomName=="BoxInBox")
	{
		this->Ph_BoxInBox=new CML2Ph_BoxInBox();
	}

// read the messenger data related to the phantom selected 
	G4UImanager* UI = G4UImanager::GetUIpointer();
	G4String command = "/control/execute ";
	UI->ApplyCommand(command+this->PhantomSpecficationsFileName); // DEVE LEGGERE TUTTI I MESSENGER DELLE MACRO STRUTTURE
}
void CML2PhantomConstruction::Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG)
{
// a call to select the right phantom
	this->design();
// create the phantom-world box
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4ThreeVector halfSize, centre;
	centre.set(0.*mm, 0.*mm, 0.*mm);
	halfSize.set(300.*mm, 300*mm, 300*mm);
	G4Box *phmWorldB = new G4Box("phmWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *phmWorldLV = new G4LogicalVolume(phmWorldB, Vacuum, "phmWorldL", 0, 0, 0);
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceWireframe(true);
	phmWorldLV->SetVisAttributes(simpleAlSVisAtt);

	G4RotationMatrix *rotation=new G4RotationMatrix();
	this->PVPhmWorld= new G4PVPlacement(rotation, centre, "phmWorldPV", phmWorldLV, PVWorld, false, 0);

// create the actual phantom
	if (this->phantomName=="fullWater")
	{
		this->Ph_fullWater->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
	}
	else if (this->phantomName=="BoxInBox")
	{
		this->Ph_BoxInBox->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
	}
}

