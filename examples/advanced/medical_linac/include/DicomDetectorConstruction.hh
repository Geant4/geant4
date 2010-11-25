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
// Author: P. Arce
// History: 30.11.07  First version
//*******************************************************
//
// DicomDetectorConstruction.hh :
//	- Start the building of the geometry
//	- Initialisation of materials
//      - Creation of the world 
//	- Reading of the DICOM data
//*******************************************************
//
// The code is part of the DICOM extended example and it was modified by :
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


#ifndef DicomDetectorConstruction_h
#define DicomDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DicomPatientZSliceHeader.hh"

 // added ML2 start
#include "G4SDManager.hh"  
#include "G4ProductionCuts.hh"  

#include "ML2SDWithParticle.hh"  
#include "ML2SDWithVoxels.hh"  
#include "ML2ReadOutGeometry.hh"
#include "ML2_DicomMessenger.hh"  

class CML2_DicomMessenger; 
// end

class G4Material;
class G4Box;
class G4LogicalVolume;

class DicomDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DicomDetectorConstruction();
  ~DicomDetectorConstruction();

  G4VPhysicalVolume *Construct(){return this->world_phys;};


  void Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG); // changed ML2
  // trigger the construction of the geometry

 // added ML2 start
  virtual void ConstructPatient()   =0;
  // construct the patient volumes. This method should be implemented for each of the derived classes
	inline CML2SDWithVoxels* getSensDet(){return  this->sensDet;};
	inline G4VPhysicalVolume *getPhysicalVolume(){return this->world_phys;};

  	inline void setDataFileName(G4String val){this->fileData=val;};
	inline void setCalibrationDensityFileName(G4String val){this->calibrationDensityFileName=val;};
	inline void setDicomDirectory(G4String val){this->dicomDirectory=val;};
	inline void setDicomColorMap(G4String val){this->dicomColorMap=val;};
	inline void setNROGVoxelsX(G4int val){this->nROGVoxelsX=val;};
	inline void setNROGVoxelsY(G4int val){this->nROGVoxelsY=val;};
	inline void setNROGVoxelsZ(G4int val){this->nROGVoxelsZ=val;};

	inline G4String getDataFileName(){return this->fileData;};
	inline G4String getCalibrationDensityFileName(){return this->calibrationDensityFileName;};
	inline G4String getDicomDirectory(){return this->dicomDirectory;};
	inline G4String getDicomColorMap(){return this->dicomColorMap;};

	inline G4int getTotalNumberOfEvents(){return this->sensDet->getTotalNumberOfEvents();};
	inline void setSurfaceToTargetZValue(G4double val){this->surfaceToTargetZValue=val;};
	inline G4String* getMatNames(){return this->matNames;};
	inline G4int getNmatNames(){return this->nMatNames;};
	G4ThreeVector getHalfContainerSize();
	void writeInfo();
private:
	G4String dicomDirectory, fileData, calibrationDensityFileName, dicomColorMap;
	CML2_DicomMessenger *dicomMessenger;

	void ML2_G4NIST_Materials();
	G4double *density; 
	G4String *inputMatNames;
	G4String *matNames;
	G4int nMatNames, nInputMatNames;
	G4ThreeVector halfContainerSize;
// end

protected:
  void InitialisationOfMaterials();
  // create the original materials

  void ReadPatientData(G4String directory, G4String fileName); // changed ML2
  // read the DICOM files describing the patient

  void ReadPatientDataFile(const G4String& fname);
  // read one of the DICOM files describing the patient (usually one per Z slice). Build a DicomPatientZSliceHeader for each file

  void MergeZSliceHeaders();
  // merge the slice headers of all the files

  G4Material* BuildMaterialWithChangingDensity( const G4Material* origMate, G4double density, G4String newMateName );
  // build a new material if the density of the voxel is different to the other voxels

  G4String ftoa(float flo);
  // convert a float to a string

  void ConstructPatientContainer();

  void SetScorer(G4LogicalVolume* voxel_logic);


protected:
  G4Material* air;

  // World ...
  G4Box* world_solid;
  G4LogicalVolume* world_logic;
  G4VPhysicalVolume* world_phys;

  G4Box* container_solid;
  G4LogicalVolume* container_logic;
  G4VPhysicalVolume* container_phys;

  G4int fNoFiles; // number of DICOM files
  std::vector<G4Material*> fOriginalMaterials;  // list of original materials 
  std::vector<G4Material*> fMaterials;  // list of new materials created to distinguish different density voxels that have the same original materials
  size_t* fMateIDs; // index of material of each voxel
 //unsigned int* fMateIDs; // index of material of each voxel

  std::map<G4int,G4double> fDensityDiffs; // Density difference to distinguish material for each original material (by index)
 
  std::vector<DicomPatientZSliceHeader*> fZSliceHeaders; // list of z slice header (one per DICOM files)
  DicomPatientZSliceHeader* fZSliceHeaderMerged; // z slice header resulted from merging all z slice headers

  G4int nVoxelX, nVoxelY, nVoxelZ;
  G4double voxelHalfDimX,  voxelHalfDimY, voxelHalfDimZ;

 // added ML2 start
CML2SDWithVoxels *sensDet;
G4int saving_in_ROG_Voxels_every_events; 
G4int seed; 
G4String ROGOutFile; 
G4bool bSaveROG;
G4int nROGVoxelsX, nROGVoxelsY, nROGVoxelsZ;
G4Region *regVol;
G4double surfaceToTargetZValue;
// end
};

#endif

