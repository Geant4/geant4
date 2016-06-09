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


#ifndef CML2ConvergenceH
#define CML2ConvergenceH

#include "G4Step.hh"

#include "ML2SinputData.hh"
#include "ML2ExpVoxels.hh"

class CML2Convergence
{
public:
	CML2Convergence();
	CML2Convergence(G4int seed, G4int saving_in_Selected_Voxels_every_events, G4String FileExperimentalData, G4String FileExperimentalDataOut, G4bool bCompareExp, G4int maxNumberOfEvents, G4int nRecycling, G4int nMaxLoops);
	~CML2Convergence(void);
	void add(const G4Step* aStep);
	G4bool stopRun();
	void setMaxNumberOfEvents(G4int val){maxNumberOfEvents=val;}
	G4double getMaxNumberOfEvents(){return maxNumberOfEvents;};
	CML2ExpVoxels * getExpVoxels(){return ML2ExpVoxels;}
	inline void saveResults(){if (ML2ExpVoxels!=0){ML2ExpVoxels->saveResults();}}
	inline void setNewGeometry(){nGeometry++;	idCurrentLoop=nMaxLoops;}
	inline int getNMaxLoops(){return nMaxLoops;}
private:
	G4bool convergenceCriteria();

	std::vector <Svoxel> voxels;
	CML2ExpVoxels *ML2ExpVoxels;

	G4String fileExperimentalData;

	G4bool bCompareExp;
        G4int maxNumberOfEvents, nGeometry, nAccumulatedEvents;
	int nMaxLoops, idCurrentLoop;
};

#endif

