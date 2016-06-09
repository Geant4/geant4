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
// ^INAIL DIPIA - ex ISPESL  and INFN Roma, gruppo collegato SanitÃ , Italy
// *Istituto Superiore di SanitÃ  and INFN Roma, gruppo collegato SanitÃ , Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2_DicomMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

CML2_DicomMessenger::CML2_DicomMessenger(DicomDetectorConstruction *dicomDetectorConstruction) : dicomDetectorConstruction(dicomDetectorConstruction)
{
	this->dicomDirectory=new G4UIcmdWithAString("/dicom/dicomDirectory",this);
	this->dicomDirectory->SetDefaultValue("");
	this->dicomDirectory->SetGuidance("directory of the dicom files");
	this->dicomDetectorConstruction->setDicomDirectory("");

	this->fileData=new G4UIcmdWithAString("/dicom/fileData",this);
	this->fileData->SetDefaultValue("");
	this->fileData->SetGuidance("file name of the dicom data.dat");
	this->dicomDetectorConstruction->setDataFileName("");

	this->calibrationDensityFileName=new G4UIcmdWithAString("/dicom/calibrationDensityFileName",this);
	this->calibrationDensityFileName->SetDefaultValue("");
	this->calibrationDensityFileName->SetGuidance("file name of the dicom CT2Density.dat");
	this->dicomDetectorConstruction->setCalibrationDensityFileName("");

	this->dicomColorMap=new G4UIcmdWithAString("/dicom/dicomColorMap",this);
	this->dicomColorMap->SetDefaultValue("");
	this->dicomColorMap->SetGuidance("file name of the dicom colormap.dat");
	this->dicomDetectorConstruction->setDicomColorMap("");

	this->nROGVoxelsX=new G4UIcmdWithAnInteger("/dicom/nROGVoxelsX",this);
	this->nROGVoxelsX->SetDefaultValue(100);
	this->nROGVoxelsX->SetGuidance("Number of voxels in the X direction");
	this->dicomDetectorConstruction->setNROGVoxelsX(100);
	
	this->nROGVoxelsY=new G4UIcmdWithAnInteger("/dicom/nROGVoxelsY",this);
	this->nROGVoxelsY->SetDefaultValue(100);
	this->nROGVoxelsY->SetGuidance("Number of voxels in the Y direction");
	this->dicomDetectorConstruction->setNROGVoxelsY(100);
	
	this->nROGVoxelsZ=new G4UIcmdWithAnInteger("/dicom/nROGVoxelsZ",this);
	this->nROGVoxelsZ->SetDefaultValue(100);
	this->nROGVoxelsZ->SetGuidance("Number of voxels in the Z direction");
	this->dicomDetectorConstruction->setNROGVoxelsZ(100);
}

CML2_DicomMessenger::~CML2_DicomMessenger(void)
{
	delete dicomDirectory;
	delete fileData;
	delete calibrationDensityFileName;
	delete dicomColorMap;
	delete nROGVoxelsX;
	delete nROGVoxelsY;
	delete nROGVoxelsZ;
}
void CML2_DicomMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==this->dicomDirectory)
	{this->dicomDetectorConstruction->setDicomDirectory(newValue);}
	else if (cmd==this->fileData)
	{this->dicomDetectorConstruction->setDataFileName(newValue);}
	else if (cmd==this->calibrationDensityFileName)
	{this->dicomDetectorConstruction->setCalibrationDensityFileName(newValue);}
	else if (cmd==this->dicomColorMap)
	{this->dicomDetectorConstruction->setDicomColorMap(newValue);}
	else if (cmd==this->nROGVoxelsX)
	{this->dicomDetectorConstruction->setNROGVoxelsX(this->nROGVoxelsX->GetNewIntValue(newValue));}
	else if (cmd==this->nROGVoxelsY)
	{this->dicomDetectorConstruction->setNROGVoxelsY(this->nROGVoxelsY->GetNewIntValue(newValue));}
	else if (cmd==this->nROGVoxelsZ)
	{this->dicomDetectorConstruction->setNROGVoxelsZ(this->nROGVoxelsZ->GetNewIntValue(newValue));}

}
