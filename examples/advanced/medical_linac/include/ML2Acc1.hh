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


#ifndef CML2Acc1H
#define CML2Acc1H

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "ML2SinputData.hh"
#include "G4ProductionCuts.hh"

class CML2Acc1Messenger;
class CML2Acc1
{
public:
	CML2Acc1(void);
	~CML2Acc1(void);
	static CML2Acc1* GetInstance(void);
	void Construct(G4VPhysicalVolume *PVWorld, G4double isoCentre);
	void reset();
	inline void setJaw1X(G4double val){this->jaw1XAperture=val;};
	inline void setJaw2X(G4double val){this->jaw2XAperture=val;};
	inline void setJaw1Y(G4double val){this->jaw1YAperture=val;};
	inline void setJaw2Y(G4double val){this->jaw2YAperture=val;};
	inline void setIsoCentre(G4double val){this->isoCentre=val;};
	inline void setidEnergy(G4int val){this->idEnergy=val;};
	inline void setLeavesAx(G4double val){this->leavesA.push_back(val);};
	inline void setLeavesBx(G4double val){this->leavesB.push_back(val);};
	inline int getidEnergy(){return this->idEnergy;};
	G4double getBeforeJaws_Z_PhaseSpacePosition(){return 215.;};
	void writeInfo();
private:
	G4double jaw1XAperture, jaw2XAperture, jaw1YAperture, jaw2YAperture, isoCentre; 
	std::vector <G4double> leavesA, leavesB;
	G4int idEnergy;
	CML2Acc1Messenger *acc1Messenger;
	static CML2Acc1 * instance;

	G4Material * otherMaterials(const G4String materialName);
	void SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4double aperture, G4RotationMatrix *cRotation);
	bool target();
	bool primaryCollimator();
	bool BeWindow();
	bool flatteningFilter();
	bool ionizationChamber();
	bool mirror();
	bool Jaw1X();
	bool Jaw2X();
	bool Jaw1Y();
	bool Jaw2Y();
	bool MLC();
	G4VPhysicalVolume * PVWorld;
};

#endif

