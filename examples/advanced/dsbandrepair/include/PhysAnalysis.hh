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
/// \file PhysAnalysis.hh
/// \brief Definition of the PhysAnalysis class

#ifndef PHYSANALYSIS_h
#define PHYSANALYSIS_h 1

#include "G4ThreeVector.hh"
#include <map>
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);               \
    void operator=(const TypeName&)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct InfoForChemGeo // store info for creating input files for chem stage
{
    G4double fType{0.}; // water: 1; solvated electron=2
    G4double fState{-99.}; // no state for solvated electron
    G4double fElectronicLevel{-99.}; // no electronic level for solvated electron
    G4double fX{0.}; // position of the incoming track
    G4double fY{0.}; // position of the incoming track
    G4double fZ{0.}; // position of the incoming track
    G4double fParentTrackID{-1.};
    G4double fEventNumber{-1.};
    G4double fVolume{-1.};
    G4double fVolumeCopyNumber{-1.};
    G4double fMotherVolume{-1.};
    G4double fMotherVolumeCopyNumber{-1.};
    G4double fRelX{-1.};
    G4double fRelY{-1.};
    G4double fRelZ{-1.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct InfoInPhysStage // store info created in Physical stage
{
    G4double fFlagParticle{-1.};
    G4double fFlagParentID{-1.};
    G4double fFlagProcess{-1.};
    G4double fX{-1.};
    G4double fY{-1.};
    G4double fZ{-1.};
    G4double fEdep{-1.};
    G4double fEventNumber{-1.};
    G4double fVolumeName{-1.};
    G4double fCopyNumber{-1.};
    G4double fLastMetVoxelCopyNum{-1.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysAnalysis
{
public:
    ~PhysAnalysis() = default;
    static PhysAnalysis* GetAnalysis();

    void OpenFile(const G4String& fname);
    void Save();
    void Close(G4bool reset = true);

    void Book();
    G4AnalysisManager* GetAnalysisManager();
    void ClearVector() ; // being called in Beginofeventaction
    void AddInfoForChemGeo(InfoForChemGeo);
    void AddInfoInPhysStage(InfoInPhysStage);
    void UpdateChemInputDataAndFillNtuple(); // being called in Endofeventaction
private:
    PhysAnalysis() = default;
    G4String CreateChemInputFile(G4int eventNum,G4int volumeCopyNumber,const G4String &voxelName);
    void UpdatingChemInputFile(InfoForChemGeo);
    void UpdatingChemInputFile(InfoInPhysStage);
    std::vector<InfoForChemGeo> fInfoForChemGeoVector;
    std::vector<InfoInPhysStage> fInfoInPhysStageVector;
    std::map<G4double, std::map<G4double, G4String> > fOutputFiles;
    G4String fOutputFolder="";
    DISALLOW_COPY_AND_ASSIGN(PhysAnalysis);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif