//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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


#ifndef CML2InputDataH
#define CML2InputDataH

#include "G4UImessenger.hh"
#include "ML2MainMessenger.hh"
#include "ML2WorldConstruction.hh"
#include "ML2PrimaryGenerationAction.hh"

class CML2MainMessenger;

class CML2CInputData
{
public:
	CML2CInputData(void);
	~CML2CInputData(void);

	inline G4bool getbOnlyVisio(){return this->bOnlyVisio;};

	inline void setbOnlyVisio(G4bool val){this->bOnlyVisio=val;};
	inline void setPhaseSpaceCentre(G4ThreeVector val){this->inputData.generalData.centrePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setPhaseSpaceHalfSize(G4ThreeVector val){this->inputData.generalData.halfSizePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setbSavePhaseSPace(G4bool val){this->inputData.generalData.bSavePhaseSpace=val;};
	inline void setbForcePhaseSpaceBeforeJaws(G4bool val){this->inputData.generalData.bForcePhaseSpaceBeforeJaws=val;};
	inline void setbStopAtPhaseSpace(G4bool val){this->inputData.generalData.bStopAtPhaseSpace=val;};
	inline void setPhaseSpaceOutFile(G4String val){this->inputData.generalData.PhaseSpaceOutFile=val;};

	inline void setbSaveROG(G4bool val){this->inputData.generalData.bSaveROG=val;};
	inline void setROGOutFile(G4String val){this->inputData.generalData.ROGOutFile=val;};

	inline void setMaxNumberOfEvents(G4int val){this->inputData.generalData.maxNumberOfEvents=val;};
	inline void setNmaxLoop(G4int val){this->inputData.generalData.nMaxLoop=val;};
	inline G4double getMaxNumberOfEvents(){return this->inputData.generalData.maxNumberOfEvents;};

	inline void setBCompareExp(G4bool val){this->inputData.generalData.bCompareExp=val;};
	inline void setFileExperimentalData(G4String val){this->inputData.generalData.fileExperimentalData=val;};
	inline void setFileExperimentalDataOut(G4String val){this->inputData.generalData.fileExperimentalDataOut=val;};

	inline void setNBeams(G4int val){this->inputData.generalData.nBeam=val;};
	inline void setNMaxParticlesInRamPlanePhaseSpace(G4int val){this->inputData.generalData.nMaxParticlesInRamPlanePhaseSpace=val;};

	inline void setSaving_in_Selected_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_Selected_Voxels_every_events=val;};
	inline void setSaving_in_ROG_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_ROG_Voxels_every_events=val;};
	inline void setMax_N_particles_in_PhSp_File(G4int val){this->inputData.generalData.max_N_particles_in_PhSp_File=val;};

	G4bool bOnlyVisio;
	SInputData inputData;
	CML2MainMessenger *ML2MainMessenger;
private:
};


#endif

