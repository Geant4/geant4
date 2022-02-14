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


#pragma once

enum DNAVolumeType
{
    deoxyribose1 = 1, 
    phosphate1, 
    deoxyribose2, 
    phosphate2,
            
    base_cytosine, 
    base_guanine, 
    base_thymine, 
    base_adenine,
            
    deoxyribose1_water, 
    phosphate1_water, 
    deoxyribose2_water, 
    phosphate2_water, 
            
    base_adenine_water, 
    base_guanine_water,
    base_cytosine_water, 
    base_thymine_water, 
            
    histone, 
    physWorld, 
    VoxelStraight
};
