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

	inline std::vector <Svoxel> getVoxels(){return this->voxels;}

	G4int getMinNumberOfEvents();
	G4int getMaxNumberOfEvents();

	G4bool loadData();

	inline void setRecycling(int recycling){this->nRecycling=recycling;};
	void saveResults(void);
private:
	void saveHeader();
	void calculateNormalizedEd(std::vector <Svoxel> &voxels);
	std::vector <Svoxel> voxels;
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

