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

#include "ML2Ph_FullWater.hh"
#include "ML2Ph_FullWaterMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4PVReplica.hh"


CML2Ph_FullWater::CML2Ph_FullWater()
{
  // phantom size and position
  halfSize.set(150.*mm,150.*mm,150.*mm);
  // phantom position
  centre.set(0.,0.,0.);
  
  fPhantomSize.setX(300.*mm);
  fPhantomSize.setY(300.*mm);
  fPhantomSize.setZ(300.*mm);
  fullWaterMessenger = new CML2Ph_FullWaterMessenger(this);
}

CML2Ph_FullWater::~CML2Ph_FullWater(void)
{
}

void CML2Ph_FullWater::writeInfo()
{
	G4cout<<"\n\n\tcentre of the phantom: " <<centre/mm<<" [mm]"<< G4endl;
	G4cout<<"\thalf thickness of the phantom: " <<halfSize/mm<<" [mm]\n"<< G4endl;
}
bool CML2Ph_FullWater::Construct(G4VPhysicalVolume *PWorld, G4int nx, G4int ny, G4int nz)
{

 PVWorld=PWorld;

 bool bCreated=false;
 G4Material *WATER=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
 G4Box *fullWaterPhantomBox = new G4Box("fullWaterPhantomBox", halfSize.getX(), halfSize.getY(), halfSize.getZ());
 G4LogicalVolume *fullWaterPhantomLV = new G4LogicalVolume(fullWaterPhantomBox, WATER, "fullWaterPhantomLV", 0, 0, 0);
 fullWaterPhantomPV = new G4PVPlacement(0, centre, "fullWaterPhantomPV", fullWaterPhantomLV, PVWorld, false, 0);

 
  G4int nxCells = nx;
  G4int nyCells = ny;
  G4int nzCells = nz;

  G4cout << "VoxelX, voxelY, Voxelz = " << nx << ", " << ny <<  ", " << nz << "; " << G4endl;

  G4ThreeVector sensSize;
  sensSize.setX(fPhantomSize.x()/(G4double)nxCells);
  sensSize.setY(fPhantomSize.y()/(G4double)nyCells);
  sensSize.setZ(fPhantomSize.z()/(G4double)nzCells);

  // The phantom is voxelised in 3D

  G4String yRepName("RepY");

  G4VSolid* solYRep = new G4Box(yRepName,fPhantomSize.x()/2.,sensSize.y()/2.,fPhantomSize.z()/2.);

  G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep,WATER,yRepName);
 
  new G4PVReplica(yRepName,logYRep,fullWaterPhantomLV,kYAxis,ny,sensSize.y());

  G4String xRepName("RepX");

  G4VSolid* solXRep = new G4Box(xRepName,sensSize.x()/2.,sensSize.y()/2.,fPhantomSize.z()/2.);
  
  G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,WATER,xRepName);
  
  new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,nx,sensSize.x());
 
  G4String zVoxName("phantomSens");
 
  G4VSolid* solVoxel = new G4Box(zVoxName,sensSize.x()/2.,sensSize.y()/2.,sensSize.z()/2.);
  
  G4LogicalVolume* LVPhantomSens = new G4LogicalVolume(solVoxel,WATER,zVoxName); // This is the Sensitive Volume

  new G4PVReplica(zVoxName,LVPhantomSens,logXRep,kZAxis,nz,sensSize.z());
//..............................................
  // Phantom segmentation using Parameterisation
  //..............................................
 
  G4cout << "  Water Phantom Size " << fPhantomSize/mm       << G4endl;
  G4cout << "  Segmentation  ("<< nx<<","<<ny<<","<<nz<<")"<< G4endl;
  
  // Region for cuts
  G4Region *regVol= new G4Region("fullWaterPhantomR");
  G4ProductionCuts* cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.1*mm);
  regVol->SetProductionCuts(cuts);

	fullWaterPhantomLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(fullWaterPhantomLV);

	// Visibility
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	fullWaterPhantomLV->SetVisAttributes(simpleAlSVisAtt);

	G4MultiFunctionalDetector* myScorer = new G4MultiFunctionalDetector("PhantomSD");
	G4SDManager::GetSDMpointer()->AddNewDetector(myScorer);
	LVPhantomSens->SetSensitiveDetector(myScorer);
       
	G4VPrimitiveScorer * totalDose = new G4PSDoseDeposit3D("TotalDose", nx,ny,nz);
	myScorer->RegisterPrimitive(totalDose);
	G4cout << "scorer registered: totalDose" << G4endl;

	bCreated=true;
	return bCreated;
}

