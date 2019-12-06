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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Updated: H. Tran (IRSN), France: 20/12/2018
//

#include "DNAParser.hh"
#include <fstream>
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Molecule struct def
struct DNAParser::Molecule
{
    Molecule(std::string name,
             int copyNumber,
             const G4ThreeVector& position,
             double radius,
             double waterRadius,
             const std::string& material,
             int strand)
                  : fName(name),
                    fMaterial(material),
                    fCopyNumber(copyNumber),
                    fPosition(position),
                    fRadius(radius),
                    fRadiusWater(waterRadius),
                    fStrand(strand)
    {}

    G4String fName;
    G4String fMaterial;
    G4int fCopyNumber;
    G4ThreeVector fPosition;
    G4double fRadius;
    G4double fRadiusWater;
    G4int fStrand;

    G4bool operator<(const Molecule& rhs) const
    {
        return (fPosition.z() < rhs.fPosition.z());
    }    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DNAParser::DNAParser()
    : fSize(0)
    , fGeoName("")
    , fpWater(nullptr)
    , fpGun(new G4MoleculeGun())   
{
    EnumParser();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAParser::~DNAParser() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DNAParser::CreateLogicVolume()
{
    G4NistManager * pMan = G4NistManager::Instance();
    fpWater = pMan->FindOrBuildMaterial("G4_WATER");

    G4String boxNameSolid = fGeoName + "_solid";
    
    G4Box* pBoxSolid = new G4Box(boxNameSolid,
                                 fSize/2,
                                 fSize/2,
                                 fSize/2);

    G4String boxNameLogic = fGeoName + "_logic";
    
    G4LogicalVolume* pBoxLogic = new G4LogicalVolume(pBoxSolid,
                                                     fpWater,
                                                     boxNameLogic);

    std::sort( fMolecules.begin(), fMolecules.end() );

    for(int i = 0, ie = fMolecules.size(); i < ie; ++i)
    {
        G4String name = fMolecules[i].fName;
        G4double radius = fMolecules[i].fRadius;
        G4double waterRadius = fMolecules[i].fRadiusWater;
        G4ThreeVector moleculePosition = fMolecules[i].fPosition;

        int copyNum = fMolecules[i].fCopyNumber;
          
        G4Orb* pMoleculeWaterSolid = nullptr;
        G4VSolid* pMoleculeWaterCutSolid = nullptr;
        G4LogicalVolume* pMoleculeWaterLogic = nullptr;
        G4VPhysicalVolume* pPhysicalVolumeWater = nullptr;   

        G4double tolerence = 1e-4 * nm;

        if(waterRadius > tolerence)
        {
            G4String nameWaterSolid = name + "_water_solid";
            pMoleculeWaterSolid = new G4Orb(nameWaterSolid, waterRadius);
            pMoleculeWaterCutSolid = CreateCutSolid(pMoleculeWaterSolid,
                                                    fMolecules[i],
                                                    fMolecules,
                                                    false);

            G4String nameWaterLogic = name + "_water_logic";
            
            pMoleculeWaterLogic = new G4LogicalVolume(pMoleculeWaterCutSolid,
                                                      fpWater,
                                                      nameWaterLogic);
            G4String nameWaterPhys = name + "_water";
            pPhysicalVolumeWater = new G4PVPlacement(0,
                                                     moleculePosition,
                                                     pMoleculeWaterLogic,
                                                     nameWaterPhys,
                                                     pBoxLogic,
                                                     false,
                                                     copyNum);
                                                     
            auto it = fEnumMap.find(nameWaterPhys);
            if(it  != fEnumMap.end())
            {
                fGeometryMap.insert(std::make_pair(pPhysicalVolumeWater,
                                                    it->second));
            }
            else
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription <<nameWaterPhys
                                     <<" could not be exsited";
                G4Exception("DNAParser::CreateLogicVolume()",
                            "DNAParser001", FatalException,
                             exceptionDescription);
            }         
        }
// Dna volume part

        G4Orb* pMoleculeSolid = nullptr;
        G4VSolid* pMoleculeCutSolid = nullptr;
        G4LogicalVolume* pMoleculeLogic = nullptr;
        G4VPhysicalVolume* pPhysicalVolumeMolecule = nullptr;   

        G4String nameSolid = fMolecules[i].fName + "_solid";
        pMoleculeSolid = new G4Orb(nameSolid, radius);
        
        pMoleculeCutSolid = CreateCutSolid(pMoleculeSolid,
                                           fMolecules[i],
                                           fMolecules,
                                           true);
       
        G4String nameLogic = name + "_logic";
        
        pMoleculeLogic = new G4LogicalVolume(pMoleculeCutSolid,
                                             fpWater,
                                             nameLogic);
        if(waterRadius > tolerence)
        {
            G4ThreeVector position(0,0,0);
            std::string namePhys = name;
            pPhysicalVolumeMolecule = new G4PVPlacement(0,
                                                        position,
                                                        pMoleculeLogic,
                                                        namePhys,
                                                        pMoleculeWaterLogic,
                                                        false,
                                                        copyNum);
                              
                              
            auto it = fEnumMap.find(namePhys);
            if(it  != fEnumMap.end())
            {
                fGeometryMap.insert(std::make_pair(pPhysicalVolumeMolecule,
                                                    it->second));
            }
            else
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription <<namePhys
                                     <<" could not be exsited";
                G4Exception("DNAParser::CreateLogicVolume()",
                            "DNAParser002", FatalException,
                             exceptionDescription);
            }                         
        }
        else
        {
            G4ThreeVector position = moleculePosition;
            G4String namePhys = name;
            pPhysicalVolumeMolecule =  new G4PVPlacement(0,
                                                         position,
                                                         pMoleculeLogic,
                                                         namePhys,
                                                         pBoxLogic,
                                                         false,
                                                         copyNum);
            
            auto it = fEnumMap.find(namePhys);
            
            if(it  != fEnumMap.end())
            {
                fGeometryMap.insert(std::make_pair(pPhysicalVolumeMolecule,
                                                    it->second));
            }
            else
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription <<namePhys
                                     <<" could not be exsited";
                G4Exception("DNAParser::CreateLogicVolume()",
                            "DNAParser003", FatalException,
                             exceptionDescription);
            }                                
        }
    }
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();
    return pBoxLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAParser::ParseFile(const std::string& fileName)
{
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();

    std::ifstream file(fileName.c_str());

    if(!file.is_open())
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription <<fileName
                             <<"could not be opened";
        G4Exception("DNAParser::ParseFile()",
                    "DNAParser002", FatalException,
                    exceptionDescription);
    }
    std::string line;
    while(std::getline(file, line) )
    {
        if( line.empty() ) continue;
        std::istringstream issLine(line);
        std::string firstItem;
        issLine >> firstItem;
        
        if(firstItem == "_Name")
        {
            G4String name = "";
            issLine >> name;
            fGeoName = name;
        }
        if(firstItem == "_Size")
        {
            G4double size;
            issLine >> size;
            size *= nm;
            fSize = size;
        }
        if(firstItem == "_Radius")
        {
            G4String name;
            issLine >> name;
            G4double radius;
            issLine >> radius;
            radius *= nm;
            G4double waterRadius;
            issLine >> waterRadius;
            waterRadius *= nm;
            fRadiusMap[name] = radius;
            fWaterRadiusMap[name] = waterRadius;
        }
        if(firstItem == "_pl")
        {
            std::string name;
            issLine >> name;
            G4String material;
            issLine >> material;

            G4int strand;
            issLine >> strand;

            G4int copyNumber;
            issLine >> copyNumber;

            G4double x;
            issLine >> x;
            x *= nm;

            G4double y;
            issLine >> y;
            y *= nm;

            G4double z;
            issLine >> z;
            z *= nm;

            Molecule molecule(name,
                              copyNumber,
                              G4ThreeVector(x, y, z),
                              fRadiusMap[name],
                              fWaterRadiusMap[name],
                              material,
                              strand);
            fMolecules.push_back(molecule);
            
            auto itt = fEnumMap.find(name);
            
            if(itt != fEnumMap.end())
            {
            
                if(itt->second != DNAVolumeType::phosphate1 && 
                   itt->second != DNAVolumeType::phosphate2)
                {
                    if(itt->second == DNAVolumeType::deoxyribose1 || 
                       itt->second == DNAVolumeType::deoxyribose2)
                    {
                        name = "Deoxyribose";
                    }
                    if(itt->second == DNAVolumeType::base_adenine)
                    {
                        name = "Adenine";
                    }
                    if(itt->second == DNAVolumeType::base_thymine)
                    {
                        name = "Thymine";
                    }
                    if(itt->second == DNAVolumeType::base_cytosine)
                    {
                        name = "Cytosine";
                    }
                    if(itt->second == DNAVolumeType::base_guanine)
                    {
                        name = "Guanine";
                    }
                    if(itt->second == DNAVolumeType::histone)
                    {
                        name = "Histone";
                    }
                
                    fpGun->AddMolecule(name, 
                                       G4ThreeVector(x, y, z), 
                                       1*picosecond);
                }
            }
            else
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription <<name
                                     <<" could not be exsited";
                G4Exception("DNAParser::ParseFile()",
                            "DNAParser005", FatalException,
                             exceptionDescription);
            }
        }
    } 
    file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSolid* DNAParser::CreateCutSolid(G4Orb *pSolidOrbRef,
                                        Molecule &molRef,
                                        std::vector<Molecule> &molList,
                                        G4bool in)
{
// The idea behing this method is to cut overlap volumes by selecting 
//one of them (the reference) and checking all the other volumes (the targets).
// If a reference and a target volumes are close enough to overlap they will be cut.
// The reference is already selected when we enter this method.
// Use the tiny space to differentiate the frontiers (may not be necessary)
    G4double tinySpace = 0.001*nm;

    G4SubtractionSolid* pSolidCut(NULL);
    G4bool isCutted = false;
    G4bool isOurVol = false;
    G4double radiusRef;

    if(molRef.fRadiusWater == 0) 
    {
        radiusRef = molRef.fRadius;
    }
    else 
    { 
        radiusRef = molRef.fRadiusWater;
    }
    
    G4ThreeVector posRef = molRef.fPosition;

    if(std::abs(posRef.x() ) + radiusRef > fSize/2 // along x
    || std::abs(posRef.y() ) + radiusRef > fSize/2 // along y
    || std::abs(posRef.z() ) + radiusRef > fSize/2) // along z
    {
        G4Box* solidBox = new G4Box("solid_box_for_cut", 
                                    fSize/2, 
                                    fSize/2, 
                                    fSize/2);
        G4ThreeVector posBox;
        G4RotationMatrix rotMat;
        rotMat.rotate(0, G4ThreeVector(0,0,1) );

        if(std::abs( posRef.y() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef +  G4ThreeVector(0,fSize,0);
            if(!isCutted)
            {
                pSolidCut = 
                new G4SubtractionSolid("solidCut",
                                       pSolidOrbRef,
                                       solidBox,
                                       &rotMat,
                                       posBox);
                isCutted = true;
            }
            else 
            {
                pSolidCut = 
                new G4SubtractionSolid("solidCut",
                                        pSolidCut,
                                        solidBox,
                                        &rotMat,
                                        posBox);
            }
        }
// Down
        if(std::abs( posRef.y() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,-fSize,0);

            if(!isCutted)
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                    pSolidOrbRef,
                                                    solidBox,
                                                    &rotMat,
                                                    posBox);
                isCutted = true;
            }
            else 
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidCut,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
            }
        }
// Left
        if(std::abs( posRef.x() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(fSize,0,0);
            if(!isCutted)
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidOrbRef,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
                isCutted = true;
            }
            else 
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidCut,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
            }
        }

// Right
        if(std::abs( posRef.x() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(-fSize,0,0);
            if(!isCutted)
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidOrbRef,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
                isCutted = true;
            }
            else 
            {   
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidCut,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
            }
        }
// Forward
        if(std::abs( posRef.z() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,fSize);
            if(!isCutted)
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidOrbRef,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
                isCutted = true;
            }
            else 
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidCut,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
            }
        }

// Backward
        if(std::abs( posRef.z() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,-fSize);
            if(!isCutted)
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidOrbRef,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
                isCutted = true;
            }
            else 
            {
                pSolidCut = new G4SubtractionSolid("solidCut",
                                                   pSolidCut,
                                                   solidBox,
                                                   &rotMat,
                                                   posBox);
            }
        }
    }

    for(int i=0, ie=molList.size(); i<ie; ++i)
    {
        G4ThreeVector posTar = molList[i].fPosition;
        G4double rTar = posRef.z();
        G4double zTar = posTar.z();

        if(zTar>rTar+20*nm)
        {
            break;
        }
        else if(zTar<rTar-20*nm)
        {
            continue;
        }

        G4double radiusTar;
        if(molList[i].fRadiusWater == 0) 
        {
            radiusTar = molList[i].fRadius;
        }
        else 
        {
            radiusTar = molList[i].fRadiusWater;
        }
     
        G4double distance = std::abs( (posRef - posTar).getR() );
      
        if(distance==0 && !isOurVol)
        {
            isOurVol = true;
            continue;
        }
        else if(distance == 0 && isOurVol)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription <<"CreateCutSolid: Two volumes "
                                 <<"are placed at the same position.";
            G4Exception("DNAParser::CreateCutSolid()",
                        "DNAParser004", FatalException,
                        exceptionDescription);
        }
        else if(distance <= radiusRef+radiusTar)
        {
            G4Box* solidBox = new G4Box("solid_box_for_cut",
                                        2*radiusTar,
                                        2*radiusTar,
                                        2*radiusTar);
            G4ThreeVector diff = posTar - posRef;
            G4double d = (std::pow(radiusRef,2) -
                          std::pow(radiusTar,2) +
                          std::pow(distance,2) )/(2*distance) +
                          solidBox->GetZHalfLength() -
                          tinySpace;
            if(in) 
            {
                d -= 2*tinySpace;
            }
            G4ThreeVector pos = d *( diff/diff.getR() );
            G4double phi = std::acos(pos.getZ()/pos.getR());
            G4double theta = std::acos( pos.getX() / 
                                ( pos.getR() * 
                                std::cos(CLHEP::pi/2. - phi) ) );
                           
            if(pos.getY()<0) 
            {
                theta = -theta;
            }
            
            G4ThreeVector rotAxisForPhi(1*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta + CLHEP::pi/2.);

            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            G4ThreeVector rotZAxis(0.,0.,1*nm);
            rotMat->rotate(theta, rotZAxis);

            if(!isCutted) 
            {
                pSolidCut = new G4SubtractionSolid("solidCut", 
                                                   pSolidOrbRef, 
                                                   solidBox, 
                                                   rotMat, 
                                                   pos);
            }
            else 
            {
            pSolidCut = new G4SubtractionSolid("solidCut", 
                                               pSolidCut, 
                                               solidBox, 
                                               rotMat, 
                                               pos);
            }
            isCutted = true;
        }
    }

    if(isCutted) 
    {
        return pSolidCut;
    }
    else 
    {
        return pSolidOrbRef;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAParser::EnumParser()
{             
    fEnumMap["deoxyribose1"]  = deoxyribose1;
    fEnumMap["phosphate1"]    = phosphate1;
    fEnumMap["deoxyribose2"]  = deoxyribose2;
    fEnumMap["phosphate2"]    = phosphate2;
    
    fEnumMap["base_cytosine"] = base_cytosine;
    fEnumMap["base_guanine"]  = base_guanine;
    fEnumMap["base_thymine"]   = base_thymine;  
    fEnumMap["base_adenine"]   = base_adenine;
    
    fEnumMap["deoxyribose1_water"] = deoxyribose1_water;
    fEnumMap["phosphate1_water"]   = phosphate1_water;
    fEnumMap["deoxyribose2_water"] = deoxyribose2_water;
    fEnumMap["phosphate2_water"]   = phosphate2_water;
    
    
    fEnumMap["base_adenine_water"]   = base_adenine_water;
    fEnumMap["base_guanine_water"]   = base_guanine_water;
    fEnumMap["base_cytosine_water"]  = base_cytosine_water;
    fEnumMap["base_thymine_water"]   = base_thymine_water;
    
    fEnumMap["histone"]  = histone;
    fEnumMap["physWorld"]    = physWorld;
    fEnumMap["VoxelStraight"] = VoxelStraight;
}

