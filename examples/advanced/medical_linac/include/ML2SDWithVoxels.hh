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


#ifndef CCML2SDWithVoxelsH
#define CCML2SDWithVoxelsH

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "ML2SinputData.hh"

class CML2ReadOutGeometryVoxels;


class CML2SDWithVoxels : public G4VSensitiveDetector
{
public:
	CML2SDWithVoxels(G4String name, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4ThreeVector centre, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ);
	~CML2SDWithVoxels(void);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist);
	inline void Initialize(G4HCofThisEvent *){};
	inline void EndOfEvent(G4HCofThisEvent*){};
	G4int getTotalNumberOfEvents(){return this->nTotalEvents;};
	inline void setActive(G4bool bActive){this->bActive=bActive;};
	void save();
	inline void setRecycling(int recycling){this->nRecycling=recycling;};
	inline void setVolumeNameIdLink(std::vector <SvolumeNameId> volumeNameIdLink){this->volumeNameIdLink=volumeNameIdLink;};
	void resetVoxelsSingle();
	void setFullOutFileDataSingle(G4String val);
private:
	void saveData(G4String Filename, Svoxel ***voxels);
	G4int getIdFromVolumeName(G4String name);

	G4ThreeVector halfSize, centre;
	G4ThreeVector pos;
	G4double halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ;
	G4int NumberOfVoxelsAlongX, NumberOfVoxelsAlongY, NumberOfVoxelsAlongZ;
	Svoxel ***voxelsSum, ***voxelsSingle;
	G4String fullOutFileData, fullOutFileDataSingle;
	G4int nTotalEvents, nSingleTotalEvents, nParticle, nParticleValatile, saving_in_ROG_Voxels_every_events, nRecycling;
	G4bool bActive, bSaveROG;
	G4double voxelMass, density, voxelVolume;
	std::vector <SvolumeNameId> volumeNameIdLink;
};



#endif
