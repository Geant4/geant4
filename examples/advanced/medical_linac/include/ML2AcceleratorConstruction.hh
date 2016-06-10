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

#ifndef CML2AcceleratorConstructionH
#define CML2AcceleratorConstructionH

#include "ML2SinputData.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4UImanager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "ML2PrimaryGenerationAction.hh"

#include "ML2Acc1.hh"

class CML2AcceleratorConstructionMessenger;
class CML2PrimaryGenerationAction;

class CML2AcceleratorConstruction
{
public:
	CML2AcceleratorConstruction(void);
	~CML2AcceleratorConstruction(void);
	static CML2AcceleratorConstruction* GetInstance(void);
	bool Construct(G4VPhysicalVolume *PVWorld, G4bool bOnlyVisio);
	inline G4VPhysicalVolume *getPhysicalVolume(void){return PVAccWorld;}
	void resetAccelerator();

	inline void setAcceleratorName(G4String val){AcceleratorName=val;}
	inline void setAcceleratorMacFileName(G4String val){AcceleratorMacFileName=val;}

	G4String getCurrentRotationString();

	inline G4String getNextAcceleratorXRotationName(){return nextAcceleratorXRotationName;}
	inline void setIsoCentre(G4double val){isoCentre=val;}
	inline void setRotation90Y(G4bool val){bRotate90Y=val;}

	inline void addAcceleratorRotationsX(G4double val){rotationsX.push_back(val);}

	inline G4double getAcceleratorIsoCentre(){return isoCentre;}
	inline G4String getAcceleratorName(){return AcceleratorName;}
	inline G4String getAcceleratorMacFileName(){return AcceleratorMacFileName;}
	inline G4double getZ_Value_PhaseSpaceBeforeJaws(){return Z_Value_PhaseSpaceBeforeJaws;}
	inline G4bool getRotation90Y(){return bRotate90Y;}
	void writeInfo();

	G4RotationMatrix * rotateAccelerator();
	G4RotationMatrix * rotateAccelerator(G4double angleX);
private:
	bool design(void);

	CML2AcceleratorConstructionMessenger *acceleratorConstructionMessenger;
	static CML2AcceleratorConstruction * instance;
	G4String AcceleratorName, AcceleratorMacFileName, nextAcceleratorXRotationName;

	G4VPhysicalVolume *PVAccWorld;
	G4int idCurrentRotationX;
	G4double currentRotationX, isoCentre, Z_Value_PhaseSpaceBeforeJaws;
	std::vector <G4double> rotationsX;
	G4ThreeVector initialCentre;
	G4bool bRotate90Y, bOnlyVisio;
	

	CML2Acc1 *accelerator1;
};

#endif
