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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
// Code based on the Geant4 extended example DICOM
//
// 

#include "ICRP110PhantomNestedParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

ICRP110PhantomNestedParameterisation::
ICRP110PhantomNestedParameterisation(const G4ThreeVector& halfVoxelSize,
                                   std::vector<G4Material*>& mat,
                                  G4int fnX_, G4int fnY_, G4int fnZ_)
:
  fdX(halfVoxelSize.x()), fdY(halfVoxelSize.y()), fdZ(halfVoxelSize.z()), 
  fnX(fnX_), fnY(fnY_), fnZ(fnZ_), // number of voxels along X, Y and Z
  fMaterials(mat), // vector of defined materials
  fMaterialIndices(nullptr) // vector which associates MaterialID to voxels
{
 ReadColourData();//Define the color of each material
                  // from ColourMap.dat
}

ICRP110PhantomNestedParameterisation::~ICRP110PhantomNestedParameterisation()
{}
 
void ICRP110PhantomNestedParameterisation::ReadColourData()
{
    // By default the tissues are not visible. Then
    // the visualisation attributes are defined based on 
    // ColourMap.dat 
    auto blankAtt = new G4VisAttributes;
    blankAtt->SetVisibility( FALSE );
    fColours["Default"] = blankAtt;

    G4String colourFile = "ColourMap.dat";
    G4cout << "Phantom Material Colours set via ColourMap.dat data file " << G4endl;

    std::ifstream fin(colourFile.c_str());
    G4int nMate;
    G4String mateName;
    G4double cred, cgreen, cblue, copacity;
    fin >> nMate;
    G4VisAttributes* visAtt=nullptr;
    
    for( G4int ii = 0; ii < nMate; ii++ ){
        fin >> mateName >> cred >> cgreen >> cblue >> copacity;
        G4Colour colour( cred, cgreen, cblue, copacity );
        visAtt = new G4VisAttributes( colour );
        visAtt->SetForceSolid(true);
        fColours[mateName] = visAtt;
       // G4cout << mateName << " colour set :  "  << colour << G4endl;
    }
}

void ICRP110PhantomNestedParameterisation::
SetNoVoxel( G4int nx, G4int ny, G4int nz )
{
  fnX = nx;
  fnY = ny;
  fnZ = nz;
}

G4Material* ICRP110PhantomNestedParameterisation::
ComputeMaterial(G4VPhysicalVolume* physVol, const G4int iz,
                const G4VTouchable* parentTouch)
{
    // protection for initialization and vis at idle state
    //
    if(parentTouch == nullptr)
        return fMaterials[0];
  
    // Copy number of voxels.
    // Copy number of X and Y are obtained from replication number.
    // Copy number of Z is the copy number of current voxel.
	
    G4int ix = parentTouch -> GetReplicaNumber(0);
    G4int iy = parentTouch -> GetReplicaNumber(1);

    G4int copyID = ix + fnX*iy + fnX*fnY*iz;
    //G4cout << "ix: "<< ix << ", iy: " << iy << ", iz:" << iz<< G4endl;
    //G4cout << "copyID from the Nested Param: "<< copyID << G4endl;
    
    //The copyID identifies the voxel    
    std::size_t matIndex = GetMaterialIndex(copyID); 
    static G4Material* mate = nullptr;
    mate = fMaterials[matIndex];
   
    if(true && physVol && G4VVisManager::GetConcreteInstance()) {
        G4String mateName = fMaterials.at(matIndex)->GetName();
        std::string::size_type iuu = mateName.find("__");
	 
        if( iuu != std::string::npos ) {
            mateName = mateName.substr( 0, iuu ); // Associate material
        }

        if(0 < fColours.count(mateName)) {
            physVol -> GetLogicalVolume() ->
            SetVisAttributes(fColours.find(mateName)->second);
        } 
        else {
            physVol->GetLogicalVolume() ->
            SetVisAttributes(fColours.begin() ->second); // Associate color
        }
    }
	physVol -> GetLogicalVolume()->SetMaterial(mate);
     
    return mate;
}

G4int ICRP110PhantomNestedParameterisation::GetMaterialIndex( G4int copyNo ) const
{
	return fMaterialIndices[copyNo];
}

G4int ICRP110PhantomNestedParameterisation::GetNumberOfMaterials() const
{
    return fMaterials.size();
}

G4Material* ICRP110PhantomNestedParameterisation::GetMaterial(G4int i) const
{
    return fMaterials[i];
}

void ICRP110PhantomNestedParameterisation::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Position of voxels.
    // x and y positions are already defined in DetectorConstruction by using
    // replicated volume. Hre we define the position along the z axis of voxels.

    physVol -> SetTranslation(G4ThreeVector(0.,0.,(2.*static_cast<double>(copyNo)
                                                +1.)*fdZ - fdZ*fnZ)); 
}

void ICRP110PhantomNestedParameterisation::
ComputeDimensions( G4Box& box, const G4int, const G4VPhysicalVolume* ) const
{
    box.SetXHalfLength(fdX);
    box.SetYHalfLength(fdY);
    box.SetZHalfLength(fdZ);
}
