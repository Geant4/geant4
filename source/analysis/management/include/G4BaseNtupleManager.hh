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

// The base class for Ntuple managers.
// It implements common functions independent from the output type.
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4BaseNtupleManager_h
#define G4BaseNtupleManager_h 1

#include "G4VNtupleManager.hh"
#include "globals.hh"

class G4BaseNtupleManager : public G4VNtupleManager
{
  public:
    explicit G4BaseNtupleManager(const G4AnalysisManagerState& state);
    G4BaseNtupleManager() = delete;
    ~G4BaseNtupleManager() override = default;

    // deleted copy constructor & assignment operator
    G4BaseNtupleManager(const G4BaseNtupleManager& rhs) = delete;
    G4BaseNtupleManager& operator=(const G4BaseNtupleManager& rhs) =delete;

  protected:
    // Methods for handling ntuples
    G4int CreateNtuple(G4NtupleBooking* booking) override = 0;

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId
    G4bool FillNtupleIColumn(G4int id, G4int value) final;
    G4bool FillNtupleFColumn(G4int id, G4float value) final;
    G4bool FillNtupleDColumn(G4int id, G4double value) final;
    G4bool FillNtupleSColumn(G4int id, const G4String& value) final;
    G4bool AddNtupleRow() final;

    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)
    G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) override = 0;
    G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) override = 0;
    G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) override = 0;
    G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId, const G4String& value) override = 0;
    G4bool AddNtupleRow(G4int ntupleId) override = 0;

    // Fisrt column id
    G4bool SetFirstNtupleColumnId(G4int firstId) final;

  protected:
    G4int   fFirstNtupleColumnId { 0 };
};

#endif

