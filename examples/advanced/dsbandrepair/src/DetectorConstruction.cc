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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "InformationKeeper.hh"
#include "DetectorConstructionMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Filesystem.hh"

RunningMode gRunMode = RunningMode::Phys;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4double factor, G4int verbose, G4bool isVisu) :
    G4VUserDetectorConstruction(), fFactor(factor), fVerbose(verbose), fBVisu(isVisu)
{
    fDetectorMessenger = new DetectorConstructionMessenger(this);
    fWorldBoxSizeX = fWorldBoxSizeY = fWorldBoxSizeZ = 1*nm;
    if (gRunMode == RunningMode::Chem) fChemGeoImport = std::make_unique<ChemGeoImport>();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()

{
    G4NistManager * man = G4NistManager::Instance();
    fWater = man->FindOrBuildMaterial("G4_WATER");
    if (gRunMode == RunningMode::Phys) ConstructFullCellNucleusGeo();
    else if (gRunMode == RunningMode::Chem) ConstructVoxelGeo();
    else {
        // only a world volume
        G4ExceptionDescription msg;
        msg <<"Only world volume is constructed."
            <<" Make sure you choose the correct running mode."
            <<" Ignore this message if you intentionally test the code.";
        G4Exception("DetectorConstruction::Construct()", 
                    "", JustWarning, msg);
        fSolidWorld = new G4Box("solidWorld",fWorldBoxSizeX, fWorldBoxSizeY, fWorldBoxSizeZ);
        fLogicWorld = new G4LogicalVolume(fSolidWorld, fWater, "logicWorld");
        fPhysWorld = new G4PVPlacement(0,G4ThreeVector(),"physWorld", 
        fLogicWorld, 0, false, false);
    }
    return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume * DetectorConstruction::ConstructFullCellNucleusGeo()
{
    PhysGeoImport geo(fBVisu);
    geo.SetFactor(fFactor);

    std::map<G4String, G4LogicalVolume*> voxelMap;
    // Create an empty logical nucleus

    G4LogicalVolume* logicNucleus = geo.CreateNucleusLogicVolume(fCellDefFilePath);
    if (!fUsingUserDefinedSizesForWorld) {
        G4double scale = 2.0;
        fWorldBoxSizeX = scale*geo.GetNucleusSizeData()["SemiX"];
        fWorldBoxSizeY = scale*geo.GetNucleusSizeData()["SemiY"];
        fWorldBoxSizeZ = geo.GetNucleusSizeData()["SemiZ"];
    }
   
    G4cout<<"===>>World box sizes: SemiX = "<<G4BestUnit(fWorldBoxSizeX,"Length")<<" , SemiY = "
          <<G4BestUnit(fWorldBoxSizeY,"Length")<<", SemiZ = "
          <<G4BestUnit(fWorldBoxSizeZ,"Length")<<"<<===="<<G4endl;
    fSolidWorld = new G4Box("solidWorld",fWorldBoxSizeX, fWorldBoxSizeY, fWorldBoxSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fWater, "logicWorld");
    fPhysWorld = new G4PVPlacement(0,G4ThreeVector(),"physWorld", fLogicWorld, 0, false, false);
    // then place cell nucleus in world
    new G4PVPlacement(0, G4ThreeVector(), logicNucleus, "nucleus_pl", fLogicWorld, false, false);
  
    G4LogicalVolume* logicStraightVoxel{nullptr}, *logicUpVoxel{nullptr}, *logicDownVoxel{nullptr}, 
                        *logicLeftVoxel{nullptr},*logicRightVoxel{nullptr};
    G4LogicalVolume* logicStraightVoxel2{nullptr},*logicUpVoxel2{nullptr},*logicDownVoxel2{nullptr}
                        , *logicLeftVoxel2{nullptr},*logicRightVoxel2{nullptr};
    
    for (const auto & entry : fVoxelDefFilesList) {
        G4String name="";
        auto ptemp = geo.CreateLogicVolume(entry,name);
        if (name=="voxelStraight" || name=="VoxelStraight") {
            logicStraightVoxel = ptemp;
            voxelMap["VoxelStraight"] = logicStraightVoxel;
        }
        if (name=="voxelUp" || name=="VoxelUp") {
            logicUpVoxel = ptemp;
            voxelMap["VoxelUp"] = logicUpVoxel;
        }
        if (name=="voxelDown" || name=="VoxelDown") {
            logicDownVoxel = ptemp;
            voxelMap["VoxelDown"] = logicDownVoxel;
        }
        if (name=="voxelRight" || name=="VoxelRight") {
            logicRightVoxel = ptemp;
            voxelMap["VoxelRight"] = logicRightVoxel;
        }
        if (name=="voxelLeft" || name=="VoxelLeft") {
            logicLeftVoxel = ptemp;
            voxelMap["VoxelLeft"] = logicLeftVoxel;
        }
        if (name=="voxelStraight2" || name=="VoxelStraight2") {
            logicStraightVoxel2 = ptemp;
            voxelMap["VoxelStraight2"] = logicStraightVoxel2;
        }
        if (name=="voxelUp2" || name=="VoxelUp2") {
            logicUpVoxel2 = ptemp;
            voxelMap["VoxelUp2"] = logicUpVoxel2;
        }
        if (name=="voxelDown2" || name=="VoxelDown2") {
            logicDownVoxel2 = ptemp;
            voxelMap["VoxelDown2"] = logicDownVoxel2;
        }
        if (name=="voxelRight2" || name=="VoxelRight2") {
            logicRightVoxel2 = ptemp;
            voxelMap["VoxelRight2"] = logicRightVoxel2;
        }
        if (name=="voxelLeft2" || name=="VoxelLeft2") {
            logicLeftVoxel2 = ptemp;
            voxelMap["VoxelLeft2"] = logicLeftVoxel2;
        }
    }

    // Create the voxel data table
    std::vector<Voxel>* voxelTable = geo.CreateVoxelsData(fCellDefFilePath);
    const G4int nucleusSize = voxelTable->size();
    G4cout<<"============================================================================"<<G4endl;
    G4cout<<"=====> Number of Histones in each voxel: "<<G4endl;
    for (auto [key, value] : geo.GetVoxelNbHistoneMap()) {
        G4cout<<key<<": \t\t\t"<<value<<G4endl;
    }
    G4cout<<"=====> Number of Basepairs in each voxel: "<<G4endl;
    for (auto [key, value] : geo.GetVoxelNbBpMap()) {
        G4cout<<key<<": \t\t\t"<<value<<G4endl;
    }
    G4cout  <<"=====> Total Number of Histones placed in geometry: \t"
            <<geo.GetTotalNbHistonePlacedInGeo()<<G4endl;
    G4cout  <<"=====> Total Number of Basepairs placed in geometry: \t"
            <<geo.GetTotalNbBpPlacedInGeo()<<G4endl;
    G4cout  <<"=====> Number of each chromatin type placed in geometry: "<<G4endl;
    for (auto [key, value] : geo.GetChromatinTypeCountMap()) {
        G4String chromatinname = "Heterochromatin";
        if (key == ChromatinType::fEuchromatin) chromatinname = "Euchromatin";
        G4cout<<chromatinname<<": \t\t\t"<<value<<G4endl;
    }
    G4cout<<"============================================================================"<<G4endl;
    
    // Create the voxel parameterisation
    if (voxelMap.size() > 0) {
        // The following dummy declarations are nescessary for SetLogicalVolum() 
        // in VoxelParameterisation to prevent from coredump
        G4double prevZpos = -fWorldBoxSizeZ;
        for (auto it = voxelMap.begin(); it != voxelMap.end(); it++) {
            G4double posX = fWorldBoxSizeX - geo.GetVoxelFullSize();
            G4double posY = fWorldBoxSizeY - geo.GetVoxelFullSize();
            G4double posZ = prevZpos + geo.GetVoxelFullSize();
            new G4PVPlacement(0, G4ThreeVector(posX,posY,posZ), 
                            it->second, "xx", fLogicWorld, false, 0);
        }
        // End dummy declaration
        G4int nthreads=0;
#ifdef G4MULTITHREADED
        nthreads = G4RunManager::GetRunManager()->GetNumberOfThreads();
#endif
        if ( (nthreads >1) && (voxelMap.size() >1)) {
            G4ExceptionDescription msg;
            msg <<"Number of thread is "<<nthreads<<" > 1; Thus, dsbandrepair will run in testing mode. "
                <<"\nThere will be no DNA constituents placed inside Cell Nucleus!!!";
            G4Exception("DetectorConstruction::ConstructFullCellNucleusGeo", 
                        "RunningMode", JustWarning, msg);
        } else {
            G4VPVParameterisation* voxelParam = new VoxelParameterisation(voxelMap, voxelTable);
            new G4PVParameterised("VoxelParam", voxelMap.begin()->second, 
                            logicNucleus, kUndefined, nucleusSize, voxelParam);
        }
    } else {
        G4ExceptionDescription msg;
        msg <<"It seems that voxel-definition files are not provided."
            <<" There will be no DNA constituents placed inside Cell Nucleus!!!";
        G4Exception("DetectorConstruction::ConstructFullCellNucleusGeo", 
                    "Geo_InputFile_NoFile", JustWarning, msg);
    }
    InformationKeeper::Instance()->RecordVoxelDefFilesList(fVoxelDefFilesList); 
    InformationKeeper::Instance()->SetTotalNbBpPlacedInGeo(geo.GetTotalNbBpPlacedInGeo());
    InformationKeeper::Instance()->SetTotalNbHistonePlacedInGeo(geo.GetTotalNbHistonePlacedInGeo());
    InformationKeeper::Instance()->SetNucleusVolume(geo.GetNucleusVolume());
    InformationKeeper::Instance()->SetNucleusMassDensity(
        logicNucleus->GetMaterial()->GetDensity()/(kg/m3)); // density in kg/m3;
    return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::ConstructVoxelGeo()
{
    if (fWorldBoxSizeX < fVoxelHalfSizeXYZ) {
        fWorldBoxSizeX = fWorldBoxSizeY = fWorldBoxSizeZ = fVoxelHalfSizeXYZ*1.1;
    }
    fSolidWorld = new G4Box("solidWorld", fWorldBoxSizeX, fWorldBoxSizeY, fWorldBoxSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fWater, "logicWorld");
    fPhysWorld = new G4PVPlacement(0,G4ThreeVector(),"physWorld", fLogicWorld, 0, false, 0);

    if (fVoxelHalfSizeXYZ>0) {
        G4cout<<"=====> Construct voxel ........"<<G4endl;
        G4Box* solidVoxel = new G4Box("solidVoxel", fVoxelHalfSizeXYZ, fVoxelHalfSizeXYZ, fVoxelHalfSizeXYZ );
        G4LogicalVolume* logicVoxel = new G4LogicalVolume(solidVoxel, fWater, "logicVoxel");
        new G4PVPlacement(0, G4ThreeVector(), logicVoxel, "physVoxel", fLogicWorld, false, 0);
    }
    return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCellDefFilePath(const G4String finput)
{
    const G4fs::path thisP{std::string(finput)};
    G4fs::file_status fst = G4fs::file_status{};
    auto isExist = G4fs::status_known(fst) ? G4fs::exists(fst) : G4fs::exists(thisP);
    if (! isExist) {
        G4String msg = "File " + finput + "does not exist !!! ";
        G4Exception("DetectorConstruction::SetNucleusDefFilePath()", 
        "Geo_InputFileNotOpened", FatalException, msg);
    } else {
        fCellDefFilePath = finput;
        InformationKeeper::Instance()->RecordCellDefFiliePath(fCellDefFilePath);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddVoxelDefFile(const G4String finput)
{
    const G4fs::path thisP{std::string(finput)};
    G4fs::file_status fst = G4fs::file_status{};
    auto isExist = G4fs::status_known(fst) ? G4fs::exists(fst) : G4fs::exists(thisP);
    if (! isExist) {
        G4String msg = "File " + finput + "does not exist !!! ";
        G4Exception("DetectorConstruction::AddVoxelDefFile()", 
        "Geo_InputFileNotOpened", FatalException, msg);
    } else {
        fVoxelDefFilesList.insert(finput);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldBoxSizes(G4ThreeVector v)
{
    G4double sizex = v.getX();
    G4double sizey = v.getY();
    G4double sizez = v.getZ();
    if (sizex > 0. && sizey > 0. && sizez > 0.) {
        fWorldBoxSizeX = sizex;
        fWorldBoxSizeY = sizey;
        fWorldBoxSizeZ = sizez;
        fUsingUserDefinedSizesForWorld = true;
    } else {
        G4ExceptionDescription msg ;
        msg << " Check your setting? One world dimensions is <= 0 !!! ";
        G4Exception("DetectorConstruction::SetWorldBoxSizes", 
        "Geo_WorldSizes", FatalException, msg);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ParseGeoFileForChemMode(const G4String fn)
{
    fChemGeoImport->ParseFiles(fn);
    fVoxelHalfSizeXYZ = fChemGeoImport->GetSize()/2.;
}