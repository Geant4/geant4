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


#ifndef CML2Ph_FullWaterH
#define CML2Ph_FullWaterH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "ML2SDWithParticle.hh"
#include "ML2SDWithVoxels.hh"
#include "ML2ReadOutGeometry.hh"

#include "G4SDManager.hh"
#include "G4ProductionCuts.hh"


class CML2Ph_FullWater
{
public:
	CML2Ph_FullWater();
	~CML2Ph_FullWater(void);
	bool Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG);
	inline G4int getTotalNumberOfEvents(){return this->sensDet->getTotalNumberOfEvents();};
	inline CML2SDWithVoxels* getSensDet(){return  this->sensDet;};
	inline G4VPhysicalVolume *getPhysicalVolume(){return this->PVWorld;};
	inline G4ThreeVector getHalfContainerSize(){return this->halfSize;};
	void writeInfo();
private:
	G4VPhysicalVolume *PVWorld;
	CML2SDWithVoxels *sensDet;
	G4ThreeVector centre, halfSize;
	G4double surfaceToTargetZValue; 
};


#endif
