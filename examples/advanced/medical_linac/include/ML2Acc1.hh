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
	inline void setJaw1X(G4double val){jaw1XAperture=val;}
	inline void setJaw2X(G4double val){jaw2XAperture=val;}
	inline void setJaw1Y(G4double val){jaw1YAperture=val;}
	inline void setJaw2Y(G4double val){jaw2YAperture=val;}
	inline void setIsoCentre(G4double val){isoCentre=val;}
	inline void setidEnergy(G4int val){idEnergy=val;}
	inline void setLeavesAx(G4double val){leavesA.push_back(val);}
	inline void setLeavesBx(G4double val){leavesB.push_back(val);}
	inline int getidEnergy(){return idEnergy;}
	G4double getBeforeJaws_Z_PhaseSpacePosition(){return 215.;}
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
        G4VPhysicalVolume *targetA_phys;
        G4VPhysicalVolume *targetB_phys;
        G4VPhysicalVolume *UpperCollimator_phys;
        G4VPhysicalVolume *CylMinusCone_phys;
        G4VPhysicalVolume *BeWTubePV;
        G4VPhysicalVolume *FFL1A_1PV;
        G4VPhysicalVolume *FFL2_1PV;
        G4VPhysicalVolume *PCUtubeW1PV;
        G4VPhysicalVolume *PCUtubeP1PV;
        G4VPhysicalVolume *PCUtubeW2PV;
        G4VPhysicalVolume *PCUtubeP2PV;
        G4VPhysicalVolume *PCUtubeW3PV;
        G4VPhysicalVolume *PCUtubeP3PV;
        G4VPhysicalVolume *MirrorTubePV;
        G4VPhysicalVolume *phVol1X;
        G4VPhysicalVolume *phVol2X;
        G4VPhysicalVolume *phVol1Y;
        G4VPhysicalVolume *phVol2Y;
        G4VPhysicalVolume *leafPhys;

};

#endif

