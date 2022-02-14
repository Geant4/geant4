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


#ifndef CML2Ph_FullWaterH
#define CML2Ph_FullWaterH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit3D.hh"
#include "ML2Ph_FullWaterMessenger.hh"

#include "G4SDManager.hh"
#include "G4ProductionCuts.hh"


class CML2Ph_FullWater
{
public:
    CML2Ph_FullWater();
    ~CML2Ph_FullWater(void);
    bool Construct(G4VPhysicalVolume *PVWorld, 
                   G4int voxelX, G4int voxelY, G4int voxelZ);
   // inline G4int getTotalNumberOfEvents(){return sensDet->getTotalNumberOfEvents();}
    inline G4VPhysicalVolume *getPhysicalVolume(){return PVWorld;}
    inline G4ThreeVector getHalfContainerSize(){return halfSize;}
    void writeInfo();

   // void SetNxVoxels(G4int val){fNx = val;}
   // void SetNyVoxels(G4int val){fNy = val;}
   // void SetNzVoxels(G4int val){fNz = val;}
   // void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) 
    //    const{ nx=fNx; ny = fNy; nz = fNz; }

private:
    CML2Ph_FullWaterMessenger *fullWaterMessenger;

    G4VPhysicalVolume *PVWorld;
    G4VPhysicalVolume *fullWaterPhantomPV;

    G4ThreeVector centre, halfSize;
    G4ThreeVector fPhantomSize;   // Size of Water Phantom
   // G4LogicalVolume* fLVPhantomSens;
};


#endif
