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


#ifndef CML2AccSaturnH
#define CML2AccSaturnH

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
#include "G4ProductionCuts.hh"

#include "ML2Accelerator.hh"
#include "ML2AccSaturnMessenger.hh"

class CML2AccSaturn : public CML2Accelerator
{
public:
 	CML2AccSaturn(void);
 	~CML2AccSaturn(void);
 	static CML2AccSaturn* GetInstance(void);

 	void Construct(G4VPhysicalVolume *PVWorld, G4double isoCentre);
	G4double getBeforeJaws_Z_PhaseSpacePosition(){return 490.;}; // dopo DOPO DOPO DEI JAWS
	void writeInfo();
    
private:
	CML2AccSaturnMessenger *accSaturnMessenger;

	static CML2AccSaturn * instance;
	void buildMaterial_SSteel();
	void buildMaterial_XC10();
	void buildMaterial_WNICU();
	void buildMaterial_Kapton();
	G4Material *getMaterial(const G4String materialName);
	G4Material *mat_XC10, *mat_WNICU, *mat_ssteel, *mat_Kapton;

	void SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4RotationMatrix *cRotation);
	bool target();
	bool primaryCollimator();
	bool vacuumWindow();
	bool flatteningFilter();
	bool ionizationChamber();
	bool Jaw1X();
	bool Jaw2X();
	bool Jaw1Y();
	bool Jaw2Y();

	G4VPhysicalVolume *PVWorld;
};

#endif

