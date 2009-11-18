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


#ifndef CML2AcceleratorConstructionH
#define CML2AcceleratorConstructionH

#include "ML2SinputData.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "ML2Acc1.hh"

class CML2Acc1;
class CML2AcceleratorConstructionMessenger;
class CML2AcceleratorConstruction
{
public:
	CML2AcceleratorConstruction(void);
	~CML2AcceleratorConstruction(void);
	static CML2AcceleratorConstruction* GetInstance(void);
	void Construct(G4VPhysicalVolume *PVWorld);
	G4VPhysicalVolume *getPhysicalVolume(void){return this->PVAccWorld;};
	inline void setAcceleratorName(G4String val){this->AcceleratorName=val;};
	inline void setAcceleratorSpecficationsFileName(G4String val){this->AcceleratorSpecficationsFileName=val;};
	inline void setAcceleratorRotationX(G4double val){this->rotationX=val;};
	inline void setAcceleratorRotationY(G4double val){this->rotationY=val;};
	inline void setAcceleratorRotationZ(G4double val){this->rotationZ=val;};
private:
	void design(void);
	CML2AcceleratorConstructionMessenger *acceleratorConstructionMessenger;
	static CML2AcceleratorConstruction * instance;
	G4String AcceleratorName, AcceleratorSpecficationsFileName;

	G4VPhysicalVolume *PVAccWorld;
	G4RotationMatrix *rotation;
	G4double rotationX, rotationY, rotationZ;
	CML2Acc1 *accelerator1;
};


#endif
