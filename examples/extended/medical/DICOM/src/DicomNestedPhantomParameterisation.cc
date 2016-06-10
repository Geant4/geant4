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
/// \file medical/DICOM/src/DicomNestedPhantomParameterisation.cc
/// \brief Implementation of the DicomNestedPhantomParameterisation class
//
// $Id: DicomNestedPhantomParameterisation.cc 76689 2013-11-14 08:43:45Z gcosmo $
//

#include "DicomNestedPhantomParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomNestedPhantomParameterisation::
DicomNestedPhantomParameterisation(const G4ThreeVector& voxelSize,
                                   std::vector<G4Material*>& mat,
                                   G4int fnZ_, G4int fnY_, G4int fnX_)
:
  //G4VNestedParameterisation(),
  fdX(voxelSize.x()), fdY(voxelSize.y()), fdZ(voxelSize.z()),
  fnX(fnX_), fnY(fnY_), fnZ(fnZ_),
  fMaterials(mat),
  fMaterialIndices(0)
{
    ReadColourData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomNestedPhantomParameterisation::~DicomNestedPhantomParameterisation()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomNestedPhantomParameterisation::ReadColourData()
{
    //----- Add a G4VisAttributes for materials not defined in file;
    /*G4VisAttributes* blankAtt = new G4VisAttributes;
    blankAtt->SetVisibility( FALSE );
    fColours["Default"] = blankAtt;

    G4String colourFile = "ColourMap.dat";

    //----- Read file
    std::ifstream fin(colourFile.c_str());
    G4int nMate;
    G4String mateName;
    G4double cred, cgreen, cblue, copacity;
    fin >> nMate;
    for( G4int ii = 0; ii < nMate; ii++ ){
        fin >> mateName >> cred >> cgreen >> cblue >> copacity;
        G4Colour colour( cred, cgreen, cblue, copacity );
        G4VisAttributes* visAtt = (copacity > 0.) ?
        (new G4VisAttributes( colour )) : 
        (new G4VisAttributes(G4VisAttributes::Invisible));
        //visAtt->SetForceSolid(true);
        fColours[mateName] = visAtt;
    }*/

    //----- Add a G4VisAttributes for materials not defined in file;
    G4VisAttributes* blankAtt = new G4VisAttributes;
    blankAtt->SetVisibility( FALSE );
    fColours["Default"] = blankAtt;

    //----- Read file
    G4String colourFile = "ColourMap.dat";
    std::ifstream fin(colourFile.c_str());
    G4int nMate;
    G4String mateName;
    G4double cred, cgreen, cblue, copacity;
    fin >> nMate;
    for( G4int ii = 0; ii < nMate; ii++ ){
        fin >> mateName >> cred >> cgreen >> cblue >> copacity;
        G4Colour colour( cred, cgreen, cblue, copacity );
        G4VisAttributes* visAtt = new G4VisAttributes( colour );
        //visAtt->SetForceSolid(true);
        fColours[mateName] = visAtt;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomNestedPhantomParameterisation::
SetNoVoxel( unsigned int nx, unsigned int ny, unsigned int nz )
{
  fnX = nx;
  fnY = ny;
  fnZ = nz;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* DicomNestedPhantomParameterisation::
ComputeMaterial(G4VPhysicalVolume* physVol, const G4int iz,
                const G4VTouchable* parentTouch)
{

    // protection for initialization and vis at idle state
    //
    if(parentTouch==0) return fMaterials[0];

    // Copy number of voxels.
    // Copy number of X and Y are obtained from replication number.
    // Copy nymber of Z is the copy number of current voxel.
    G4int ix = parentTouch->GetReplicaNumber(0);
    G4int iy = parentTouch->GetReplicaNumber(1);

    G4int copyID = ix + fnX*iy + fnX*fnY*iz;

    unsigned int matIndex = GetMaterialIndex(copyID);
    static G4Material* mate = 0;
    mate = fMaterials[matIndex];


    if(false && physVol && G4VVisManager::GetConcreteInstance()) {
        G4String mateName = fMaterials.at(matIndex)->GetName();
        std::string::size_type iuu = mateName.find("__");
        if( iuu != std::string::npos ) {
            mateName = mateName.substr( 0, iuu );
        }

        if(0 < fColours.count(mateName)) {
            physVol->GetLogicalVolume()->
            SetVisAttributes(fColours.find(mateName)->second);
        } else {
            physVol->GetLogicalVolume()->
            SetVisAttributes(fColours.begin()->second);
        }
    }
    
    return mate;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int DicomNestedPhantomParameterisation::
GetMaterialIndex( unsigned int copyNo ) const
{
    //return *(fMaterialIndices+copyNo);
    return fMaterialIndices[copyNo];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Number of Materials
// Material scanner is required for preparing physics tables and so on before
// starting simulation, so that G4 has to know number of materials.
//
G4int DicomNestedPhantomParameterisation::GetNumberOfMaterials() const
{
    return fMaterials.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// GetMaterial
//  This is needed for material scanner and realizing geometry.
//
G4Material* DicomNestedPhantomParameterisation::GetMaterial(G4int i) const
{
    return fMaterials[i];
}

//
// Transformation of voxels.
//
void DicomNestedPhantomParameterisation::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Position of voxels.
    // x and y positions are already defined in DetectorConstruction by using
    // replicated volume. Here only we need to define is z positions of voxels.
    physVol->SetTranslation(G4ThreeVector(0.,0.,(2.*static_cast<double>(copyNo)
                                                +1.)*fdZ - fdZ*fnZ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Dimensions are always same in this RE02 example.
//
void DicomNestedPhantomParameterisation::
ComputeDimensions( G4Box& box, const G4int, const G4VPhysicalVolume* ) const
{
    box.SetXHalfLength(fdX);
    box.SetYHalfLength(fdY);
    box.SetZHalfLength(fdZ);
}
