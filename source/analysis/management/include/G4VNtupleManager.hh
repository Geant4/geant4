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

// The pure abstract base class for Ntuple manager.
// It defines functions independent from the output type.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VNtupleManager_h
#define G4VNtupleManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "globals.hh"

#include <vector>

struct G4NtupleBooking;

class G4VNtupleManager : public G4BaseAnalysisManager
{
  // Disable using the object managers outside G4VAnalysisManager and
  // its messenger
  friend class G4VAnalysisManager;

  public:
    explicit G4VNtupleManager(const G4AnalysisManagerState& state)
      : G4BaseAnalysisManager(state) {}
    G4VNtupleManager() = delete;
    ~G4VNtupleManager() override = default;

    // deleted copy constructor & assignment operator
    G4VNtupleManager(const G4VNtupleManager& rhs) = delete;
    G4VNtupleManager& operator=(const G4VNtupleManager& rhs) = delete;

  protected:
    // Methods for handling ntuples
    virtual G4int CreateNtuple(G4NtupleBooking* booking) = 0;

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId
    virtual G4bool FillNtupleIColumn(G4int id, G4int value) = 0;
    virtual G4bool FillNtupleFColumn(G4int id, G4float value) = 0;
    virtual G4bool FillNtupleDColumn(G4int id, G4double value) = 0;
    virtual G4bool FillNtupleSColumn(G4int id, const G4String& value) = 0;
    virtual G4bool AddNtupleRow() = 0;

    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)
    virtual G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) = 0;
    virtual G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) = 0;
    virtual G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) = 0;
    virtual G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId,
                                     const G4String& value) = 0;
    virtual G4bool AddNtupleRow(G4int ntupleId) = 0;

    // Activation option
    virtual void  SetActivation(G4bool activation) = 0;
    virtual void  SetActivation(G4int id, G4bool activation) = 0;
    virtual G4bool  GetActivation(G4int id) const = 0;

    // New cycle option
    virtual void SetNewCycle(G4bool value) = 0;
    virtual G4bool GetNewCycle() const = 0;

    // Access methods
    virtual G4int GetNofNtuples() const = 0;
    // virtual G4int GetNofNtupleBookings() const = 0;

    // Set first column Id
    virtual G4bool SetFirstNtupleColumnId(G4int firstId) = 0;

    // Clear all date
    virtual void Clear() = 0;

    // List ntuples
    virtual G4bool List(std::ostream& output, G4bool onlyIfActive = true) = 0;
};

#endif
