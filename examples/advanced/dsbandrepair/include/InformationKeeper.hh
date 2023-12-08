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
/// \file InformationKeeper.hh
/// \brief Definition of the InformationKeeper class

#ifndef InformationKeeper_h
#define InformationKeeper_h 1

#include "globals.hh"
#include <map>
#include <set>
class InformationKeeper
{
public:
    static InformationKeeper* Instance();
    ~InformationKeeper() = default;
    void RecordPhysDNAName(const G4String &p) {fPhysDNAName = p;};
    void RecordCellDefFiliePath(const G4String &pth) {fCellDefFilePath = pth;};
    void RecordVoxelDefFilesList(std::set<G4String> list) {fVoxelDefFilesList = list;};
    void RecordChemInputFolderName(const G4String &pth) {fChemInputFolderName = pth;};
    void WritePhysGeo();
    G4String GetChemInputFolderName() {return fChemInputFolderName;}
    G4String GetPhysOutFolderName() {return fPhysOutFolderName;}
    void SetTotalNbBpPlacedInGeo(unsigned long long val) {fTotalNbBpPlacedInGeo = val;}
    void SetTotalNbHistonePlacedInGeo(unsigned long long val) {fTotalNbHistonePlacedInGeo = val;}
    void SetNucleusVolume(G4double vl) {fNucleusVolume = vl;};
    void SetNucleusMassDensity(G4double md) {fNucleusMassDensity = md;};
private:
    explicit InformationKeeper();
    static InformationKeeper* fInstance;
    
    G4String fPhysDNAName{""};
    G4String fCellDefFilePath;
    std::set<G4String> fVoxelDefFilesList;
    G4String fChemInputFolderName{""};
    G4String fPhysOutFolderName{""};
    unsigned long long fTotalNbBpPlacedInGeo{0};
    unsigned long long fTotalNbHistonePlacedInGeo{0};
    G4double fNucleusVolume{0.};
    G4double fNucleusMassDensity{0.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif