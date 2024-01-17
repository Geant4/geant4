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
/// \file ChemNtupleManager.hh
/// \brief Definition of the ChemNtupleManager class

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChemNtupleManager
{
public:
    ChemNtupleManager() = default;
    ~ChemNtupleManager() = default;

    void SetFileName(const G4String& name) {fFileName = name;};
    void Book();
    void Save();
    void FillNtupleIColumn(G4int icol, G4int ival);
    void FillNtupleFColumn(G4int icol, G4float ival);
    void FillNtupleDColumn(G4int id, G4int icol, G4double ival);
    void AddNtupleRow();
    void AddNtupleRow(G4int id);
    void CreateNtuples();

private:
    G4String fFileName="Output";// output
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
