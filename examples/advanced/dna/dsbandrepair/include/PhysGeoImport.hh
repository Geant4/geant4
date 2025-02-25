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
/// \file PhysGeoImport.hh
/// \brief Definition of the PhysGeoImport class

#ifndef GEOIMPORT_HH
#define GEOIMPORT_HH

#include <map>
#include <fstream>
#include <algorithm>
#include <array>
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4IntersectionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct Molecule
{
    Molecule(G4String name, G4int copyNumber, G4ThreeVector position, 
    G4double radius, G4double waterRadius, G4String material, G4int strand)
    {
        fName = name;
        fMaterial = material;
        fCopyNumber = copyNumber;
        fPosition = position;
        fRadius = radius;
        fRadiusWater = waterRadius;
        fStrand = strand;
    }

    G4String fName{""};
    G4String fMaterial{""};

    G4int fCopyNumber{-1};
    G4int fStrand{-1};

    G4ThreeVector fPosition;

    G4double fRadius{0.};
    G4double fRadiusWater{0.};

    // To sort the molecules in function of their z coordinate
    G4bool operator<(const Molecule& str) const
    {
        return (fPosition.z() < str.fPosition.z() );
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct Voxel
{
    enum VoxelType
    {
	Straight,
        Left,
        Right,
        Up,
        Down,
        Straight2,
        Left2,
        Right2,
        Up2,
        Down2,
        Other
    };

    Voxel(G4int copyNumber, G4int chromoNum, G4int domainNum, 
    VoxelType type, const G4ThreeVector& pos, G4RotationMatrix* rot)
    {
        fCopyNumber = copyNumber;
        fChromoNum = chromoNum;
        fDomainNum = domainNum;
        fType = type;
        fPos = pos;
        fpRot = rot;
    }

    G4int fCopyNumber{0};
    G4int fChromoNum{0};
    G4int fDomainNum{0};
    G4ThreeVector fPos;
    G4RotationMatrix* fpRot{nullptr};
    VoxelType fType{VoxelType::Other};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

enum ChromatinType{
    fUnspecified = 0,
    fHeterochromatin = 1,
    fEuchromatin = 2,
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysGeoImport
{
public:
    PhysGeoImport();
    PhysGeoImport(G4bool isVisu);
    ~PhysGeoImport() = default;

    void SetFactor(G4double factor){fFactor=factor;}

    G4double GetFactor() const {return fFactor;}

    G4String GetGeoName() const {return fGeoName;}

    // This method will trigger the parse of the file and the build of the geometry
    G4LogicalVolume* CreateLogicVolume(const G4String& fileName, G4String& voxelName);

    G4LogicalVolume* CreateNucleusLogicVolume(const G4String& fileName);

    std::vector<Voxel>* CreateVoxelsData(const G4String& fileName);
    std::map<G4String, G4int> GetVoxelNbHistoneMap() {return fVoxelNbHistoneMap;}
    std::map<G4String, G4int> GetVoxelNbBpMap() {return fVoxelNbBpMap;}
    unsigned long long GetTotalNbBpPlacedInGeo() {return fTotalNbBpPlacedInGeo;}
    unsigned long long GetTotalNbHistonePlacedInGeo() {return fTotalNbHistonePlacedInGeo;}
    G4double GetNucleusVolume() {return fNucleusVolume;}
    G4double GetVoxelFullSize() {return fSize;}
    std::map<G4String, G4double> GetNucleusSizeData() {return fNucleusData;}
    std::map<ChromatinType, unsigned long long> GetChromatinTypeCountMap() {return fChromatinTypeCount;}
private:

    G4bool fIsVisu{false};

    // Factor to scale the geometry
    G4double fFactor{1.};

    G4double fSize{0.};
    G4double fNucleusVolume{0.};
    unsigned long long fTotalNbBpPlacedInGeo = 0;
    unsigned long long fTotalNbHistonePlacedInGeo = 0;
    G4String fGeoName="";

    G4String fNucleusName="CellNucleus";
    G4String fNucleusType="";
    std::map<G4String, G4double> fNucleusData;

    std::map<G4String, G4double> fRadiusMap;
    std::map<G4String, G4double> fWaterRadiusMap;
    std::map<G4String, G4int> fVoxelNbBpMap;
    std::map<G4String, G4int> fVoxelNbHistoneMap;

    // Vector to contain all the molecule structures listed within the imput file
    std::vector<Molecule> fMolecules;

    // To check if this is the first voxel of the chromosome
    std::map<G4int, G4bool> fFirstMap;

    // Materials
    std::vector<G4Material*> fMaterialVect;
    G4Material* fpWater{nullptr};
    G4Material* fTHF{nullptr};
    G4Material* fPY{nullptr};
    G4Material* fPU{nullptr};
    G4Material* fTMP{nullptr};
    G4Material* fSugarMixt{nullptr};
    G4Material* fDeoxyribose{nullptr};
    G4Material* fPhosphate{nullptr};
    G4Material* fCytosine_PY{nullptr};
    G4Material* fThymine_PY{nullptr};
    G4Material* fGuanine_PU{nullptr};
    G4Material* fAdenine_PU{nullptr};
    G4Material* fHomogeneous_dna{nullptr};
    G4Material* fVacuum{nullptr};
    G4String  ParseFile(const G4String& fileName);
    G4VSolid* CreateCutSolid(G4Orb *solidOrbRef,
                             Molecule &molRef,
                             std::vector<Molecule> &molList,
                             G4bool in);
    void DefineMaterial();

    std::map<ChromatinType, unsigned long long> fChromatinTypeCount;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // GEOIMPORT_HH
