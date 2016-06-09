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


#include "ML2PhantomConstruction.hh"
#include "ML2PhantomConstructionMessenger.hh"

CML2PhantomConstruction::CML2PhantomConstruction(void): PVPhmWorld(0), sensDet(0)
{
	this->phantomContstructionMessenger=new CML2PhantomConstructionMessenger(this);
	this->idCurrentCentre=0;
}

CML2PhantomConstruction::~CML2PhantomConstruction(void)
{
	if (this->phantomName=="fullWater")
	{
		delete Ph_fullWater;
	}
	else if (this->phantomName=="boxInBox")
	{
		delete Ph_BoxInBox;
	}
	else  if (this->phantomName=="Dicom1")
	{
		delete Ph_Dicom;
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
bool CML2PhantomConstruction::design(void)
{
// switch between two different phantoms according to the main macro selection
	bool bPhanExists=false;

	std::cout << "I'm building "<< this->phantomName<<"  phantom"<< G4endl;

	if (this->phantomName=="fullWater")
	{
		this->Ph_fullWater=new CML2Ph_FullWater();bPhanExists=true;
		this->halfPhantomInsideSize=this->Ph_fullWater->getHalfContainerSize();
	}
	else if (this->phantomName=="boxInBox")
	{
		this->Ph_BoxInBox=new CML2Ph_BoxInBox();bPhanExists=true;
		this->halfPhantomInsideSize=this->Ph_BoxInBox->getHalfContainerSize();
	}
	else  if (this->phantomName=="Dicom1")
	{
		this->Ph_Dicom=new RegularDicomDetectorConstruction();bPhanExists=true;

	// read the messenger data related to the phantom selected 
		G4UImanager* UI = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		UI->ApplyCommand(command+this->PhantomFileName ); 

		DicomHandler *dcmHandler=new DicomHandler();
		dcmHandler->CheckFileFormat(this->Ph_Dicom->getDicomDirectory(), this->Ph_Dicom->getDataFileName(), this->Ph_Dicom->getCalibrationDensityFileName());

		this->halfPhantomInsideSize=this->Ph_Dicom->getHalfContainerSize();
	}

 	if (this->centre.size()<1)
 	{this->addNewCentre(G4ThreeVector(0.,0.,0.));}
	return bPhanExists;
}
G4int CML2PhantomConstruction::getTotalNumberOfEvents()
{
	if (this->phantomName="fullWater")
	{return this->Ph_fullWater->getTotalNumberOfEvents();}
	else if (this->phantomName="boxInBox")
	{return this->Ph_BoxInBox->getTotalNumberOfEvents();}
	else if (this->phantomName="Dicom1")
	{return this->Ph_Dicom->getTotalNumberOfEvents();}
	return 0;
}

bool CML2PhantomConstruction::Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4bool bOnlyVisio)
{
	this->idVolumeName=0;
	this->bOnlyVisio=bOnlyVisio;
// a call to select the right phantom
	if(this->design())
	{
		this->phantomContstructionMessenger->SetReferenceWorld(bOnlyVisio);	// create the phantom-world box
		G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

		G4Box *phmWorldB = new G4Box("phmWorldG", this->halfPhantomInsideSize.getX(), this->halfPhantomInsideSize.getY(), this->halfPhantomInsideSize.getZ());
		G4LogicalVolume *phmWorldLV = new G4LogicalVolume(phmWorldB, Vacuum, "phmWorldL", 0, 0, 0);
		G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
		simpleAlSVisAtt->SetVisibility(false);
// 		simpleAlSVisAtt->SetForceWireframe(false);
		phmWorldLV->SetVisAttributes(simpleAlSVisAtt);
	

		this->PVPhmWorld= new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "phmWorldPV", phmWorldLV, PVWorld, false, 0);
	
	// create the actual phantom
		if (this->phantomName=="fullWater")
		{
			this->Ph_fullWater->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
			this->sensDet=this->Ph_fullWater->getSensDet();
			this->createPhysicalVolumeNamesList(this->Ph_fullWater->getPhysicalVolume());
			this->Ph_fullWater->writeInfo();
		}
		else if (this->phantomName=="boxInBox")
		{
			this->Ph_BoxInBox->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
			this->sensDet=this->Ph_BoxInBox->getSensDet();
			this->createPhysicalVolumeNamesList(this->Ph_BoxInBox->getPhysicalVolume());
			this->Ph_BoxInBox->writeInfo();
		}
		else  if (this->phantomName=="Dicom1")
		{
			std::cout << this->Ph_Dicom->getDicomDirectory()<<" "<< this->Ph_Dicom->getDataFileName()<<" "<< this->Ph_Dicom->getCalibrationDensityFileName()<<" "<< this->Ph_Dicom->getDicomColorMap() << G4endl;

			this->Ph_Dicom->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
			this->sensDet=this->Ph_Dicom->getSensDet();
			this->createPhysicalVolumeNamesList(this->Ph_Dicom->getMatNames(), this->Ph_Dicom->getNmatNames());
			this->Ph_Dicom->writeInfo();
		}
		// I create the data base volumeName-volumeID in the sensitive detector 

		this->sensDet->setVolumeNameIdLink(this->volumeNameIdLink);
	}
	else
	{
		return false;
	}
	return true;
}
void CML2PhantomConstruction::createPhysicalVolumeNamesList(G4String  *matNames, G4int nMatNames)
{
	SvolumeNameId svnid;
	for (int i=0;i< nMatNames; i++)
	{
		svnid.volumeId=i;
		svnid.volumeName=matNames[i];
		this->volumeNameIdLink.push_back(svnid);
	}
}
void CML2PhantomConstruction::createPhysicalVolumeNamesList(G4VPhysicalVolume  *PV)
{
	int nLVD1;
	nLVD1=(int) PV->GetLogicalVolume()->GetNoDaughters();
	SvolumeNameId svnid;
		std::cout << "PV in name: " <<PV->GetName() << G4endl;
	if (nLVD1>0)
	{
		for (int i=0; i <nLVD1; i++)
		{
			this->createPhysicalVolumeNamesList(PV->GetLogicalVolume()->GetDaughter(i));
		}
		this->idVolumeName++;
		svnid.volumeId=this->idVolumeName;
		svnid.volumeName=PV->GetLogicalVolume()->GetMaterial()->GetName();
		this->volumeNameIdLink.push_back(svnid);
		std::cout << "physical volume name: " <<svnid.volumeName << G4endl;
	}
	else
	{
		this->idVolumeName++;
		svnid.volumeId=this->idVolumeName;
		svnid.volumeName=PV->GetLogicalVolume()->GetMaterial()->GetName();
		this->volumeNameIdLink.push_back(svnid);
		std::cout << "physical volume name: " <<svnid.volumeName << G4endl;
	}
}
bool  CML2PhantomConstruction::applyNewCentre()
{
	if (this->idCurrentCentre <(int) this->centre.size())
	{
		this->currentCentre=this->centre[this->idCurrentCentre];
		this->applyNewCentre(this->currentCentre);
		this->idCurrentCentre++;
		return true;
	}
	return false;
}
void CML2PhantomConstruction::writeInfo()
{
	if (!this->bOnlyVisio)
	{std::cout <<"Actual centre: "<<this->idCurrentCentre<<"/"<<this->centre.size() <<"  "<< G4endl;}
	std::cout <<"Phantom and its ROG centre: " << this->currentCentre<< G4endl;
}
void CML2PhantomConstruction::applyNewCentre(G4ThreeVector centre)
{
	if (this->sensDet!=0)
	{	
		this->currentCentre=centre;
		G4GeometryManager::GetInstance()->OpenGeometry();
		this->PVPhmWorld->SetTranslation(centre);
		this->sensDet->GetROgeometry()->GetROWorld()->GetLogicalVolume()->GetDaughter(0)->SetTranslation(centre);
		this->sensDet->resetVoxelsSingle();
		G4GeometryManager::GetInstance()->CloseGeometry();
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
}
G4String CML2PhantomConstruction::getCurrentTranslationString()
{
	char cT[5];
	G4int cTI;
	G4String translationName;
	cTI=(G4int)((this->currentCentre.getX()/mm));
	sprintf(cT,"%d",cTI);
	translationName="_TrX"+G4String(cT)+"_";
	cTI=(G4int)((this->currentCentre.getY()/mm));
	sprintf(cT,"%d",cTI);
	translationName+="Y"+G4String(cT)+"_";
	cTI=(G4int)((this->currentCentre.getZ()/mm));
	sprintf(cT,"%d",cTI);
	translationName+="Z"+G4String(cT);
	return translationName;
}

