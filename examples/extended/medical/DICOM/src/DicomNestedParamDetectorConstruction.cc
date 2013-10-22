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
// $Id$
//
/// \file medical/DICOM/src/DicomNestedParamDetectorConstruction.cc
/// \brief Implementation of the DicomNestedParamDetectorConstruction class
//
// History:
//        Pedro Arce
//
//*******************************************************

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "DicomNestedParamDetectorConstruction.hh"
#include "DicomNestedPhantomParameterisation.hh"

#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomNestedParamDetectorConstruction::DicomNestedParamDetectorConstruction()
 : DicomDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomNestedParamDetectorConstruction::~DicomNestedParamDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomNestedParamDetectorConstruction::ConstructPhantom()
{
#ifdef G4VERBOSE
    G4cout << "DicomNestedParamDetectorConstruction::ConstructPhantom " << G4endl;
#endif

    /*
    //----- Replication of Water Phantom Volume.
    //--- Z Slice
    G4String zRepName("RepZ");
    G4VSolid* solZRep = new G4Box(zRepName, fNVoxelX*fVoxelHalfDimX, fNVoxelY*fVoxelHalfDimY, 
    fVoxelHalfDimZ);
    G4LogicalVolume* logZRep = new G4LogicalVolume(solZRep, fAir, zRepName);
    new G4PVReplica(zRepName, logZRep, fContainer_logic, kZAxis, fNVoxelZ, fVoxelHalfDimZ*2.);

    logZRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));

    //--- X Slice
    G4String xRepName("RepX");
    G4VSolid* solXRep = new G4Box(xRepName, fNVoxelX*fVoxelHalfDimX, fVoxelHalfDimY, 
    fVoxelHalfDimZ);
    G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep, fAir, xRepName);
    new G4PVReplica(xRepName, logXRep, logZRep, kYAxis, fNVoxelY, fVoxelHalfDimY*2.);

    logXRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));

    //----- Voxel solid and logical volumes
    //--- Y Slice
    G4VSolid* solVoxel = new G4Box("phantom",fVoxelHalfDimX,fVoxelHalfDimY,fVoxelHalfDimZ);
    G4LogicalVolume* logicVoxel = new G4LogicalVolume(solVoxel,fAir,"phantom");

    //logicVoxel->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));

    //
    // Parameterisation for transformation of voxels.
    //  (voxel size is fixed in this example.
    //    e.g. nested parameterisation handles material and transfomation of voxels.)
    G4ThreeVector voxelSize(fVoxelHalfDimX,fVoxelHalfDimY,fVoxelHalfDimZ);
    DicomNestedPhantomParameterisation* param = new DicomNestedPhantomParameterisation(voxelSize,
    fMaterials);

    new G4PVParameterised("phantom",
                          logicVoxel,
                          logXRep,
                          kXAxis,
                          fNVoxelX,
                          param);*/


    //----- Replication of Water Phantom Volume.
    //--- Y Slice
    G4String yRepName("RepY");
    G4VSolid* solYRep = new G4Box(yRepName,fNVoxelX*fVoxelHalfDimX,fVoxelHalfDimY,
                  fNVoxelZ*fVoxelHalfDimZ);
    G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep,fAir,yRepName);
    new G4PVReplica(yRepName,logYRep,fContainer_logic,kYAxis,fNVoxelY,fVoxelHalfDimY*2.);

    logYRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));

    //--- X Slice
    G4String xRepName("RepX");
    G4VSolid* solXRep = new G4Box(xRepName,fVoxelHalfDimX,fVoxelHalfDimY,fNVoxelZ*fVoxelHalfDimZ);
    G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,fAir,xRepName);
    new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,fNVoxelX,fVoxelHalfDimX*2.);

    logXRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));
    
    //----- Voxel solid and logical volumes
    //--- Z Slice
    G4VSolid* solVoxel = new G4Box("phantom",fVoxelHalfDimX,fVoxelHalfDimY,fVoxelHalfDimZ);
    G4LogicalVolume* logicVoxel = new G4LogicalVolume(solVoxel,fAir,"phantom");

    logicVoxel->SetVisAttributes(new G4VisAttributes(G4VisAttributes::Invisible));

    //
    // Parameterisation for transformation of voxels.
    //  (voxel size is fixed in this example.
    //    e.g. nested parameterisation handles material and transfomation of voxels.)
    G4ThreeVector voxelSize(fVoxelHalfDimX,fVoxelHalfDimY,fVoxelHalfDimZ);
    DicomNestedPhantomParameterisation* param = new DicomNestedPhantomParameterisation(voxelSize,
                                               fMaterials);

    new G4PVParameterised("phantom",    // their name
                          logicVoxel, // their logical volume
                          logXRep,      // Mother logical volume
                          kZAxis,       // Are placed along this axis
                          //kUndefined,
                          // Are placed along this axis
                          fNVoxelZ,      // Number of cells
                          param);       // Parameterisation.


     

    param->SetMaterialIndices( fMateIDs );
    param->SetNoVoxel( fNVoxelX, fNVoxelY, fNVoxelZ );

    //phantom_phys->SetRegularStructureId(0);

    // Z logical volume
    SetScorer(logicVoxel);

    // X logical volume
    //SetScorer(logXRep);

    // Y logical volume
    //SetScorer(logYRep);

    // Container logical volume
    //SetScorer(fContainer_logic);

}


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4SDManager.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSDoseDeposit3D.hh"

void DicomNestedParamDetectorConstruction::ConstructSDandField()
{

    G4cout << "\n\n\n\n\t CONSTRUCT SD AND FIELD \n\n\n" << G4endl;

    //G4SDManager* SDman = G4SDManager::GetSDMpointer();

    //SDman->SetVerboseLevel(1);

    //
    // Sensitive Detector Name
    G4String concreteSDname = "phantomSD";
    std::vector<G4String> scorer_names;
    scorer_names.push_back(concreteSDname);
    //------------------------
    // MultiFunctionalDetector
    //------------------------
    //
    // Define MultiFunctionalDetector with name.
    // declare MFDet as a MultiFunctionalDetector scorer
    G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(concreteSDname);
    //SDman->AddNewDetector( MFDet );                 // Register SD to SDManager
    G4VPrimitiveScorer* dosedep = new G4PSDoseDeposit3D("DoseDeposit", fNVoxelX, fNVoxelY,
    fNVoxelZ);
    MFDet->RegisterPrimitive(dosedep);

    for(std::set<G4LogicalVolume*>::iterator ite = scorers.begin(); ite != scorers.end(); ++ite) {
        SetSensitiveDetector(*ite, MFDet);
    }

 //if(DicomRunAction::Instance()->GetDicomRun()) {
 //  DicomRunAction::Instance()->GetDicomRun()->ConstructMFD(scorer_names);
 //  }
    
    
}*/


