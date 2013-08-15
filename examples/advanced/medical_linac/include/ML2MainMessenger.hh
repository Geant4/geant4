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


#ifndef CML2MainMessengerH
#define CML2MainMessengerH

#include "ML2Main.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "ML2CInputData.hh"

class CML2Main;
class CML2CInputData;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAInteger;

class CML2MainMessenger : public G4UImessenger 
{
public:
	CML2MainMessenger(CML2CInputData *CInputData);
	~CML2MainMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	CML2CInputData *CInputData;

	G4UIcmdWith3VectorAndUnit *phaseSpaceCentre, *phaseSpaceHalfSize;
	G4UIcmdWithAString *phaseSPaceOutFile, *ROGOutFile;
	G4UIcmdWithABool *bSavePhaseSpace, *bStopAtPhaseSpace, *bSaveROG, *bForcePhaseSpaceBeforeJaws;

	G4UIcmdWithAnInteger *nBeam, *nMaxParticlesInRamPlanePhaseSpace, *maxNumberOfEvents;
	G4UIcmdWithABool *bCompareExp, *bOnlyVisio;
	G4UIcmdWithAString * fileExperimentalData, * fileExperimentalDataOut;
	G4UIcmdWithAnInteger *saving_in_Selected_Voxels_every_events;
	G4UIcmdWithAnInteger *saving_in_ROG_Voxels_every_events;
	G4UIcmdWithAnInteger *max_N_particles_in_PhSp_File, *nMaxLoop; 
};

#endif
