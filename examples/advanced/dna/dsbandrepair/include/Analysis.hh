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
/// \file Analysis.hh
/// \brief Definition of the Analysis class
/// \file Analysis.hh
/// \brief Definition of the Analysis class

#ifndef ANALYSIS_h
#define ANALYSIS_h 1

#include "G4ThreeVector.hh"
#include <map>
#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);               \
    void operator=(const TypeName&)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct InfoForChemGeo // store info for creating input files for chem stage
{
    G4int fType{0}; // water: 1; solvated electron=2
    G4int fState{-99}; // no state for solvated electron
    G4int fElectronicLevel{-99}; // no electronic level for solvated electron
    G4double fX{0.}; // position of the incoming track
    G4double fY{0.}; // position of the incoming track
    G4double fZ{0.}; // position of the incoming track
    G4int fParentTrackID{-1};
    G4int fEventNumber{-1};
    G4int fVolume{-1};
    G4int fVolumeCopyNumber{-1};
    G4int fMotherVolume{-1};
    G4int fMotherVolumeCopyNumber{-1};
    G4double fRelX{-1.};
    G4double fRelY{-1.};
    G4double fRelZ{-1.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct InfoInPhysStage // store info created in Physical stage
{
    G4int fFlagParticle{-1};
    G4int fFlagParentID{-1};
    G4int fFlagProcess{-1};
    G4double fX{-1.};
    G4double fY{-1.};
    G4double fZ{-1.};
    G4double fEdep{-1.};
    G4int fEventNumber{-1};
    G4int fVolumeName{-1};
    G4int fCopyNumber{-1};
    G4int fLastMetVoxelCopyNum{-1};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Analysis
{
public:
    ~Analysis() = default;
    static Analysis* GetAnalysis();

    void OpenFile(const G4String outFolder="");
    void Save();
    void Close(G4bool reset = true);
    void SetFileName(const G4String& name) {fFileName = name;};
    void Book();
    G4AnalysisManager* GetAnalysisManager();
    void ClearVector() ; // being called in Beginofeventaction
    void AddInfoForChemGeo(InfoForChemGeo);
    void AddInfoInPhysStage(InfoInPhysStage);
    void UpdateChemInputDataAndFillNtuple(); // being called in Endofeventaction
    void RecordCellDefFiliePath(const G4String &pth) {fCellDefFilePath = pth;};
    void RecordVoxelDefFilesList(std::set<G4String> list) {fVoxelDefFilesList = list;};
    void RecordChemInputFolderName(const G4String &pth) {fChemInputFolderName = pth;};
    void WritePhysGeo();
    G4String GetChemInputFolderName() {return fChemInputFolderName;}
    G4String GetPhysOutFolderName() {return fPhysOutFolderName;}
    G4String GetChemOutFolderName() {return fChemOutFolderName;}
    void SetTotalNbBpPlacedInGeo(unsigned long long val) {fTotalNbBpPlacedInGeo = val;}
    void SetTotalNbHistonePlacedInGeo(unsigned long long val) {fTotalNbHistonePlacedInGeo = val;}
    void SetNucleusVolume(G4double vl) {fNucleusVolume = vl;};
    void SetNucleusMassDensity(G4double md) {fNucleusMassDensity = md;};
    void CheckAndCreateNewFolderInChemStage();
    void CheckAndCreateNewFolderInPhysStage();
private:
    Analysis() {DefineCommands();};
    G4String CreateChemInputFile(G4int eventNum,G4int volumeCopyNumber,const G4String &voxelName);
    void UpdatingChemInputFile(InfoForChemGeo);
    void UpdatingChemInputFile(InfoInPhysStage);
    std::vector<InfoForChemGeo> fInfoForChemGeoVector;
    std::vector<InfoInPhysStage> fInfoInPhysStageVector;
    std::map<G4double, std::map<G4double, G4String> > fOutputFiles;
    G4String fCellDefFilePath;
    std::set<G4String> fVoxelDefFilesList;
    G4String fChemInputFolderName{"chem_input"};
    G4String fPhysOutFolderName{"phys_output"};
    G4String fChemOutFolderName{"chem_output"};
    unsigned long long fTotalNbBpPlacedInGeo{0};
    unsigned long long fTotalNbHistonePlacedInGeo{0};
    G4double fNucleusVolume{0.};
    G4double fNucleusMassDensity{0.};
    G4String fFileName="Output";// output
    std::unique_ptr<G4GenericMessenger> fMessenger;
    void DefineCommands();
    DISALLOW_COPY_AND_ASSIGN(Analysis);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif