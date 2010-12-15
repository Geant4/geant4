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


#include "ML2WorldConstruction.hh"

CML2WorldConstruction::CML2WorldConstruction():acceleratorEnv(0),phantomEnv(0),PVWorld(0),phaseSpace(0),backScatteredPlane(0)
{
	this->phantomEnv=CML2PhantomConstruction::GetInstance();
	this->acceleratorEnv=CML2AcceleratorConstruction::GetInstance();
	this->bWorldCreated=false;
}

CML2WorldConstruction::~CML2WorldConstruction(void)
{
	delete this->PVWorld;
	delete this->phantomEnv;
	delete this->acceleratorEnv;
	delete this->phaseSpace;
	delete this->backScatteredPlane;
}

CML2WorldConstruction* CML2WorldConstruction::instance = 0;

CML2WorldConstruction* CML2WorldConstruction::GetInstance()
{
	if (instance == 0)
	{
		instance = new CML2WorldConstruction();
	}
	return instance;
}

G4VPhysicalVolume* CML2WorldConstruction::Construct()
{
	return this->PVWorld;
}
bool CML2WorldConstruction::create(SInputData *inputData, bool bOnlyVisio)
{
// create the world box 
	this->bOnlyVisio=bOnlyVisio;
	G4double halfSize=3000.*mm;
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *worldB = new G4Box("worldG", halfSize, halfSize, halfSize);
	G4LogicalVolume *worldLV = new G4LogicalVolume(worldB, Vacuum, "worldL", 0, 0, 0);
	G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour::Black());
	simpleWorldVisAtt->SetVisibility(false);
// 	simpleWorldVisAtt->SetForceSolid(false);
	worldLV->SetVisAttributes(simpleWorldVisAtt);
	this->PVWorld= new G4PVPlacement(0,  G4ThreeVector(0.,0.,0.), "worldPV", worldLV, 0, false, 0);

// create the accelerator-world box 
	if (!this->acceleratorEnv->Construct(this->PVWorld, bOnlyVisio))
	{
		std::cout <<"\n\n The macro file '"<<inputData->generalData.StartFileInputData<<"' refers to a not defined accelerator.\n"<< this->acceleratorEnv->getAcceleratorName()<<"\n\nSTOP\n\n" << G4endl;
		return false;
	}

// create the phantom-world box 
	if (!this->phantomEnv->Construct(this->PVWorld, inputData->generalData.saving_in_ROG_Voxels_every_events, inputData->generalData.seed, inputData->generalData.ROGOutFile, inputData->generalData.bSaveROG, bOnlyVisio))
	{
		std::cout <<"\n\n The macro file '"<<inputData->generalData.StartFileInputData<<"' refers to a not defined phantom.\n"<< this->phantomEnv->getPhantomName()<<"\n\nSTOP\n\n" << G4endl;
		return false;
	}

// if the bSavePhaseSpace flag is true create a phase plane 
	if (inputData->generalData.bSavePhaseSpace)
	{
		this->phaseSpace=new CML2PhaseSpaces();
		if (inputData->generalData.bForcePhaseSpaceBeforeJaws)
		{inputData->generalData.centrePhaseSpace.setZ(this->acceleratorEnv->getZ_Value_PhaseSpaceBeforeJaws());}
		this->phaseSpace->createPlane(idSD_PhaseSpace, inputData->generalData.max_N_particles_in_PhSp_File, inputData->generalData.seed, inputData->generalData.nMaxParticlesInRamPlanePhaseSpace, this->acceleratorEnv->getPhysicalVolume(), "PhSp", inputData->generalData.PhaseSpaceOutFile, inputData->generalData.bSavePhaseSpace, inputData->generalData.bStopAtPhaseSpace, inputData->generalData.centrePhaseSpace, inputData->generalData.halfSizePhaseSpace,&inputData->primaryParticleData, this->acceleratorEnv->getAcceleratorIsoCentre());
	}

// create a killer plane to destroy the particles back scattered from the target
	this->backScatteredPlane=new CML2PhaseSpaces();
	this->backScatteredPlane->createPlane(this->acceleratorEnv->getPhysicalVolume(), "killerPlane", G4ThreeVector(0, 0, -50*mm), G4ThreeVector(200*mm, 200*mm, 1*mm));

	this->bWorldCreated=true;
	return true;
}
void CML2WorldConstruction::checkVolumeOverlap()
{
// loop inside all the daughters volumes
	std::cout<< G4endl;
	bool bCheckOverlap=false;
	int nWorld;
	nWorld=(int) this->PVWorld->GetLogicalVolume()->GetNoDaughters();

	int nSubWorlds, nSubWorlds2;
	for (int i=0; i<(int) this->PVWorld->GetLogicalVolume()->GetNoDaughters(); i++)
	{
		bCheckOverlap=this->PVWorld->GetLogicalVolume()->GetDaughter(i)->CheckOverlaps();
		nSubWorlds=(int) this->PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetNoDaughters();
		for (int j=0; j<nSubWorlds; j++)
		{
			bCheckOverlap=this->PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->CheckOverlaps();
			nSubWorlds2=(int) this->PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetNoDaughters();
			for (int k=0; k<nSubWorlds2; k++)
			{
				bCheckOverlap=this->PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetDaughter(k)->CheckOverlaps();
			}
		}
	}
	std::cout<< G4endl;
}
bool CML2WorldConstruction::newGeometry()
{
	G4bool bNewRotation=false;
	G4bool bNewCentre=false;
	G4bool bNewGeometry=false;
	bNewCentre=this->phantomEnv->applyNewCentre();
	G4RotationMatrix *rmInv=this->acceleratorEnv->rotateAccelerator();
	if (rmInv!=0)
	{
		CML2PrimaryGenerationAction::GetInstance()->setRotation(rmInv);
		bNewRotation=true;
	}
	if (bNewRotation || bNewCentre){bNewGeometry=true;}
	return bNewGeometry;
}

