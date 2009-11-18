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
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
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
	CML2Convergence(G4int seed, G4int saving_in_Selected_Voxels_every_events, G4String FileExperimentalData, G4bool bCompareExp, G4int minNumberOfEvents);
	~CML2Convergence(void);
	void add(const G4Step* aStep);
	G4bool runAgain();
private:
	G4bool convergenceCriteria();

	std::vector <Svoxel> voxels;
	CML2ExpVoxels *ML2ExpVoxels;

	G4String fileExperimentalData;

	G4bool bCompareExp;
	G4int minNumberOfEvents;
};

#endif

