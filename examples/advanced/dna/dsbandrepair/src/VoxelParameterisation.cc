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
//
/// \file VoxelParameterisation.cc
/// \brief Implementation of the VoxelParameterisation class

#include "VoxelParameterisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelParameterisation::VoxelParameterisation(std::map<G4String, 
                        G4LogicalVolume*>& voxelMap, std::vector<Voxel>* voxels)
    : G4VPVParameterisation(),
    fVoxelMap(voxelMap),
    fVoxels(voxels)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelParameterisation::~VoxelParameterisation()
{
    if(fVoxels!=0)
        delete fVoxels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Retrive the correct voxel data for the cpN
    Voxel& voxelData = fVoxels->at(copyNo);

    // Set the current physical volume parameters
    physVol->SetTranslation(voxelData.fPos);
    physVol->SetRotation(voxelData.fpRot);
    physVol->SetName(VoxelName(voxelData.fType) );
    physVol->SetLogicalVolume(LogicalVoxel(voxelData.fType) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* VoxelParameterisation::LogicalVoxel(Voxel::VoxelType type) const
{
    G4LogicalVolume* voxel = 0;

    if(type==Voxel::Straight)
        voxel = fVoxelMap.at("VoxelStraight");
    else if(type==Voxel::Right)
        voxel =fVoxelMap.at("VoxelRight");
    else if(type==Voxel::Left)
        voxel =fVoxelMap.at("VoxelLeft");
    else if(type==Voxel::Up)
        voxel =fVoxelMap.at("VoxelUp");
    else if(type==Voxel::Down)
        voxel =fVoxelMap.at("VoxelDown");
    else if(type==Voxel::Straight2)
        voxel = fVoxelMap.at("VoxelStraight2");
    else if(type==Voxel::Right2)
        voxel =fVoxelMap.at("VoxelRight2");
    else if(type==Voxel::Left2)
        voxel =fVoxelMap.at("VoxelLeft2");
    else if(type==Voxel::Up2)
        voxel =fVoxelMap.at("VoxelUp2");
    else if(type==Voxel::Down2)
        voxel =fVoxelMap.at("VoxelDown2");
    else
    {
        G4ExceptionDescription msg;
        msg << "Voxel type "<<type<<" is not registered";
        G4Exception("VoxelParameterisation::GetLogicalVoxel", "", FatalException, msg);
    }

    if(voxel==0)
    {
        G4ExceptionDescription msg;
        msg << "Voxel is a nullptr";
        G4Exception("VoxelParameterisation::GetLogicalVoxel", "", FatalException, msg);
    }

    return voxel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String VoxelParameterisation::VoxelName(Voxel::VoxelType type) const
{
    G4String name ("");

    if(type==Voxel::Straight)
        name = "VoxelStraight";
    else if(type==Voxel::Right)
        name = "VoxelRight";
    else if(type==Voxel::Left)
        name = "VoxelLeft";
    else if(type==Voxel::Up)
        name = "VoxelUp";
    else if(type==Voxel::Down)
        name = "VoxelDown";
    else if(type==Voxel::Straight2)
        name = "VoxelStraight2";
    else if(type==Voxel::Right2)
        name = "VoxelRight2";
    else if(type==Voxel::Left2)
        name = "VoxelLeft2";
    else if(type==Voxel::Up2)
        name = "VoxelUp2";
    else if(type==Voxel::Down2)
        name = "VoxelDown2";
    else
    {
        G4ExceptionDescription msg;
        msg << "Voxel type "<<type<<" is not registered";
        G4Exception("VoxelParameterisation::GetLogicalVoxel", "", FatalException, msg);
    }

    return name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


