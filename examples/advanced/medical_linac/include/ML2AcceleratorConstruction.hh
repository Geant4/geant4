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
