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
/// \file ChemGeoImporthh
/// \brief Definition of the ChemGeoImport class

#ifndef ChemGeoImport_HH
#define ChemGeoImport_HH

#include <map>
#include <fstream>
#include <algorithm>
#include <set>

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4H2O.hh"
#include "G4Electron_aq.hh"
#include "G4Scheduler.hh"

#include "UserMoleculeGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ChemMolecule
{
    ChemMolecule(G4String name, G4int copyNumber, G4ThreeVector position,
                G4int strand, G4int state, G4int electronicLevel, G4int trackId)
    {
        fName = name;
        fCopyNumber = copyNumber;
        fPosition = position;
        fStrand = strand;
        fState = state;
        fElectronicLevel = electronicLevel;
        fTrackId = trackId;
    }

    ~ChemMolecule() {}

    G4String fName{""};

    G4int fCopyNumber{-1};
    G4int fStrand{-1};
    G4int fState{-99};
    G4int fElectronicLevel{-99};
    G4int fTrackId{-1};

    G4ThreeVector fPosition{0};

    friend G4bool operator==(const ChemMolecule& lhs, const ChemMolecule& rhs)
    {
        return (lhs.fName == rhs.fName
                && lhs.fCopyNumber == rhs.fCopyNumber
                && lhs.fStrand == rhs.fStrand
                && lhs.fState == rhs.fState
                && lhs.fElectronicLevel == rhs.fElectronicLevel
                && lhs.fTrackId == rhs.fTrackId);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChemGeoImport
{
public:
    ChemGeoImport();
    ~ChemGeoImport();

    void SetFactor(double factor){fFactor=factor;}
    G4double GetFactor() const {return fFactor;}

    G4double GetSize() const {return fSize;}

    void ParseFiles(const G4String& chemInputFile);

    // This method will trigger the build of the geometry
    void InsertMoleculeInWorld();

    void Reset();
    G4String GetVoxelDefFilePath(G4String bareName);
    G4bool IsFileParsed() {return fIsParsed;}
private:
    G4bool fIsParsed{false};

    // Factor to scale the geometry
    G4double fFactor{1};

    G4double fSize{0};

    G4String fGeoNameFromChemInput{""};

    // Vector to contain all the molecule structures listed within the imput file
    std::vector<ChemMolecule> fMolecules;

    std::vector<ChemMolecule> fToBeRemovedMol;

    UserMoleculeGun* fpGun{nullptr};

    void ParseChemInputFile(const G4String& fileName);
    void ParseGeoFile(const G4String& fileName);
    G4bool IsMoleculeInTheRemoveTable(const ChemMolecule& molecule);

    void GetVoxelDefFilePathList();
    std::set<G4String> fVoxelDefFilesList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // ChemGeoImport_HH
