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
/// @file Analysis.hh
/// @brief Define histograms

#ifndef ANALYSIS_MANAGER_H
#define ANALYSIS_MANAGER_H

#include "G4ThreeVector.hh"
#include <tools/histo/h1d>
#include <tools/histo/h2d>


#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class Analysis {
public:
  ~Analysis();

  static Analysis* GetAnalysis();

  void Book();
  void EndOfRun();

  void OpenFile(const G4String& fname);
  void Save();
  void Close(G4bool reset = true);

  void FillIncident(const G4ThreeVector& p);
  void FillDose(const G4ThreeVector& p, G4double dedx);

  void ClearIncidentFlag();

  void SetUseNtuple(G4bool useNtuple) { 
    G4cout << "Set useNtuple: " << useNtuple << G4endl;
    fUseNtuple = useNtuple; 
  }

  void SetMergeNtuple(G4bool mergeNtuple) { 
    G4cout << "Set mergeNtuple: " << mergeNtuple << G4endl;
    fMergeNtuple = mergeNtuple; 
  }

private:
  Analysis();
  DISALLOW_COPY_AND_ASSIGN(Analysis);

  G4bool fUseNtuple;
  G4bool fMergeNtuple;

  //Histograms handlers
  G4int fincident_x_hist;
  G4int fincident_map;
  G4int fdose_hist;
  G4int fdose_map;
  G4int fdose_prof;
  G4int fdose_map_prof;
  G4int fdose_map3d;

  static G4ThreadLocal G4int fincidentFlag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void Analysis::ClearIncidentFlag()
{
  fincidentFlag = false;
}

#endif
