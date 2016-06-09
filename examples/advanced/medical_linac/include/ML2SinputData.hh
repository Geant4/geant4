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

#ifndef inputDataH
#define inputDataH


#include "globals.hh"
#include <vector>
#include "G4RotationMatrix.hh"

enum idParticleSource
{
	id_randomTarget=1,
	id_phaseSpace=2
};
enum idTypeOfSensitiveDetector
{
	idSD_ComponentROG=1,
	idSD_PhaseSpace=2,
	idSD_KillerPlane=3
};
struct SStartInputData
{
	G4String fileInputData;
	G4int seed;
};
struct SGeneralData
{
	G4String WorldName, fileExperimentalData, fileExperimentalDataOut, StartFileInputData;
	G4int seed, nBeam, nMaxParticlesInRamPlanePhaseSpace;
	G4bool bSaveROG, bCompareExp;
	G4String PhaseSpaceOutFile, ROGOutFile;
	G4bool bForcePhaseSpaceBeforeJaws;
	G4bool bSavePhaseSpace;
	G4bool bStopAtPhaseSpace;
	G4ThreeVector centrePhaseSpace, halfSizePhaseSpace;

	G4int maxNumberOfEvents, nMaxLoop;
	int saving_in_Selected_Voxels_every_events;
	int saving_in_ROG_Voxels_every_events;
	int max_N_particles_in_PhSp_File; 
};

struct Sparticle
{
	G4ThreeVector pos, dir;
	G4double kinEnergy;
	G4int nPrimaryPart, primaryParticlePDGE, partPDGE, volumeId;
	G4String volumeName;
};
struct SPrimaryParticle
{
	G4int partPDGE, nPrimaryParticle;
};
struct SInputData
{
	SGeneralData generalData;
	SPrimaryParticle primaryParticleData;
};
struct Svoxel
{
	G4ThreeVector pos, halfSize;
	G4double depEnergy, depEnergy2, expDose, depEnergyNorm, depEnergyNormError;
	G4int nEvents, volumeId;
};
struct SvolumeNameId
{
	G4String volumeName;
	G4int volumeId;
};
#endif
