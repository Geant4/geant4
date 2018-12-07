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

#ifndef CML2AcceleratorH
#define CML2AcceleratorH

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

class CML2Accelerator
{
public:
	CML2Accelerator(void);
	virtual ~CML2Accelerator(void);
 	virtual void Construct(G4VPhysicalVolume *PVWorld, G4double iso);
	void reset();
    void setJaw1X(G4double val);
    void setJaw2X(G4double val);
    void setJaw1Y(G4double val);
    void setJaw2Y(G4double val);
    void setIsoCentre(G4double val);
    void setidEnergy(G4int val);
    void setLeavesAx(G4double val);
    void setLeavesBx(G4double val);
    int getidEnergy();
	virtual G4double getBeforeJaws_Z_PhaseSpacePosition(){return 0;}
	virtual void writeInfo(){return;};

	G4double jaw1XAperture, jaw2XAperture, jaw1YAperture, jaw2YAperture;
	G4double isoCentre;
	std::vector <G4double> vec_leavesA, vec_leavesB;
	G4int idEnergy;
protected:


private:

};

#endif

