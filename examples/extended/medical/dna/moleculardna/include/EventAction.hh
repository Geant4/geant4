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
/// file:EventAction.hh
/// brief:
#ifndef MOLECULAR_EVENT_ACTION_HH
#define MOLECULAR_EVENT_ACTION_HH

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"
#include <vector>
#include <map>

class AnalysisManager;

class DNAHit;

class ChromosomeHit;

class DNAGeometry;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
 public:
  explicit EventAction(AnalysisManager*);

  ~EventAction() override;

  void BeginOfEventAction(const G4Event*) override;

  void EndOfEventAction(const G4Event*) override;

  void AddDNAHit(const DNAHit* hit) { fDNAHits.push_back(hit); };

  void AddChromosomeEdep(const G4String&, const G4double&, const G4bool&);

  void AddCellEdep(const G4double&);

  void AddTrackLengthCell(const G4double&);

  void AddTrackLengthChro(const G4double&);

  void SetPrimStopPos(const G4ThreeVector& pos) { fprimstoppos = pos; };

  void Initialize();

 private:
  G4bool fInitialized = false;
  AnalysisManager* fAnalysisManager;
  DNAGeometry* fDNAGeometry = nullptr;
  std::map<uint32_t, ChromosomeHit*> fChromoHitMap;
  std::vector<const DNAHit*> fDNAHits;
  G4double fEdepCell         = 0;  // dousatsu
  G4double fTraLenCell       = 0;  // dousatsu
  G4double fTraLenChro       = 0;  // dousatsu
  G4ThreeVector fprimstoppos = G4ThreeVector(0., 0., 0.);
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_EVENT_ACTION_HH
