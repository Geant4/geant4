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


#ifndef CML2ReadOutGeometryVoxelsH
#define CML2ReadOutGeometryVoxelsH

#include "G4ThreeVector.hh"
#include "ML2SinputData.hh"
#include "G4Step.hh"

class CML2ExpVoxels 
{
public:
	CML2ExpVoxels(G4bool bHasExperimentalData, G4int saving_in_Selected_Voxels_every_events, G4int seed, G4String FileExperimentalData, G4String FileExperimentalDataOut);
	~CML2ExpVoxels(void);
	void add(G4ThreeVector pos, G4double depEnergy, G4double density);
	void add(const G4Step* aStep);

	inline std::vector <Svoxel> getVoxels(){return vec_voxels;}
	inline void setRecycling(int recycling){nRecycling = recycling;}
	void saveResults(void);
	void resetNEventsInVoxels();

	G4int getMinNumberOfEvents();
	G4int getMaxNumberOfEvents();
	G4bool loadData();


private:
	void saveHeader();
	void calculateNormalizedEd(std::vector <Svoxel> &vec_voxels);

	std::vector <Svoxel> vec_voxels;
	G4int *nVoxelsgeometry;
	G4ThreeVector minZone, maxZone;
	G4int nCurves;
	G4int *startCurve, *stopCurve;
	G4double *chi2Factor;
	G4String headerText1, headerText2, fullFileIn, fullFileOut;
	G4String seedName, loopName;
	SGeneralData *generalData;
	G4int nParticle;
	G4int nTotalEvents, saving_in_Selected_Voxels_every_events, nRecycling;
	G4bool bHasExperimentalData;

};

#endif

