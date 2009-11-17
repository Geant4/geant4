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
#include "ML2SinputData.h"
#include "G4ProductionCuts.hh"

class CML2Acc1Messenger;
class CML2Acc1
{
public:
	CML2Acc1(void);
	~CML2Acc1(void);
	static CML2Acc1* GetInstance(void);
	void Construct(G4VPhysicalVolume *PVWorld);
	inline void setJaw1X(G4double val){this->jaw1XAperture=val;};
	inline void setJaw2X(G4double val){this->jaw2XAperture=val;};
	inline void setJaw1Y(G4double val){this->jaw1YAperture=val;};
	inline void setJaw2Y(G4double val){this->jaw2YAperture=val;};
	inline void setSSD(G4double val){this->SSD=val;};
	inline void setidEnergy(G4int val){this->idEnergy=val;};
	inline void setLeavesAx(G4double val){this->leavesA.push_back(val);};
	inline void setLeavesBx(G4double val){this->leavesB.push_back(val);};

private:
	G4double jaw1XAperture, jaw2XAperture, jaw1YAperture, jaw2YAperture, SSD; 
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
	G4VPhysicalVolume *PVWorld;
};

#endif

