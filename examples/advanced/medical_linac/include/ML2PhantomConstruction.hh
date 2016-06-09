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


#ifndef CML2PhantomConstructionH
#define CML2PhantomConstructionH

#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "ML2SinputData.hh"
#include "ML2Ph_FullWater.hh"
#include "ML2Ph_BoxInBox.hh"



// DICOM
#include "DicomHandler.hh"
#include "RegularDicomDetectorConstruction.hh"

class CML2PhantomConstructionMessenger;


class CML2PhantomConstruction
{
public:
	CML2PhantomConstruction(void);
	~CML2PhantomConstruction(void);
	static CML2PhantomConstruction* GetInstance(void);
	bool Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4bool bOnlyVisio);
	G4int getTotalNumberOfEvents();
	inline G4String getPhantomName(){return this->phantomName;};
	inline void setPhantomName(G4String val){this->phantomName=val;};
	inline void setPhantomFileName (G4String val){this->PhantomFileName =val;};
	inline void setNewName(){this->sensDet->setFullOutFileDataSingle("");};
	inline void setNewName(G4String val){this->sensDet->setFullOutFileDataSingle(val);};

	void applyNewCentre(G4ThreeVector val); 
	bool applyNewCentre(); // it opens the geometry changes the phantom centre and closes the geometry : used also by CML2PhantomConstructionMessenger

	inline void saveData(){this->sensDet->save();};

	inline void addNewCentre(G4ThreeVector val){this->centre.push_back(val);};

	void writeInfo();
	G4String getCurrentTranslationString();
	inline void resetSensDet(){this->sensDet->resetVoxelsSingle();};
private:
	bool design(void);
	void createPhysicalVolumeNamesList(G4String  *matNames, G4int nMatNames);
	void createPhysicalVolumeNamesList(G4VPhysicalVolume  *PV);
	CML2PhantomConstructionMessenger *phantomContstructionMessenger;
	static CML2PhantomConstruction * instance;
	G4int nVoxelsX, nVoxelsY, nVoxelsZ;
	G4String phantomName, PhantomFileName ;

	G4VPhysicalVolume *PVPhmWorld;


	CML2SDWithVoxels *sensDet;
	std::vector <SvolumeNameId> volumeNameIdLink;
	G4int idVolumeName;

	G4ThreeVector halfPhantomInsideSize, currentCentre;
	std::vector <G4ThreeVector> centre;
	G4int idCurrentCentre;

	CML2Ph_FullWater *Ph_fullWater;
	CML2Ph_BoxInBox *Ph_BoxInBox;
	DicomDetectorConstruction  *Ph_Dicom;
	G4bool bOnlyVisio;
};
#endif

