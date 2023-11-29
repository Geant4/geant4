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

// Class template for ntuple managers for all output types.
//
// Author: Ivana Hrivnacova, 19/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4TNtupleManager_h
#define G4TNtupleManager_h 1

#include "G4BaseNtupleManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include <vector>
#include <string_view>

// NT - ntuple type, FT - file type
template <typename NT, typename FT>
class G4TNtupleManager : public G4BaseNtupleManager {

  public:
    explicit G4TNtupleManager(const G4AnalysisManagerState& state);
    G4TNtupleManager() = delete;
    ~G4TNtupleManager() override;

  protected:
    // Methods to manipulate ntuples
    G4int CreateNtuple(G4NtupleBooking* ntupleBooking) override;

    virtual void CreateNtuplesFromBooking(
                   const std::vector<G4NtupleBooking*>& ntupleBookings);

    virtual G4bool Reset();
    void Clear() override;

    // Methods to create ntuples
    // are implemented in G4NtupleBookingManager base class

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId (from base class)
    using G4BaseNtupleManager::FillNtupleIColumn;
    using G4BaseNtupleManager::FillNtupleFColumn;
    using G4BaseNtupleManager::FillNtupleDColumn;
    using G4BaseNtupleManager::FillNtupleSColumn;
    using G4BaseNtupleManager::AddNtupleRow;
    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)
    G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) final;
    G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) final;
    G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) final;
    G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId, const G4String& value) final;
    G4bool AddNtupleRow(G4int ntupleId) override;

    // Activation option
    //
    void SetActivation(G4bool activation) final;
    void SetActivation(G4int ntupleId, G4bool activation) final;
    G4bool GetActivation(G4int ntupleId) const final;

    // New cycle option
    void SetNewCycle(G4bool value) override;
    G4bool GetNewCycle() const override;

    // Access methods
    NT* GetNtuple() const;
    NT* GetNtuple(G4int ntupleId) const;
    G4int GetNofNtuples() const final;

    // List ntuples
    G4bool List(std::ostream& output, G4bool onlyIfActive = true) final;

    // Iterators
    typename std::vector<NT*>::iterator BeginNtuple();
    typename std::vector<NT*>::iterator EndNtuple();
    typename std::vector<NT*>::const_iterator BeginConstNtuple() const;
    typename std::vector<NT*>::const_iterator EndConstNtuple() const;

    // Data members
    std::vector<G4TNtupleDescription<NT, FT>*> fNtupleDescriptionVector;
    std::vector<NT*> fNtupleVector;
    const std::vector<G4NtupleBooking*>* fNtupleBookingVector { nullptr };
    G4bool fNewCycle { false };

  private:
    // Methods

    // Fuctions which are specific to output type
    //
    virtual void CreateTNtupleFromBooking(
                    G4TNtupleDescription<NT, FT>* ntupleDescription) = 0;

    virtual void FinishTNtuple(
                    G4TNtupleDescription<NT, FT>* ntupleDescription,
                    G4bool fromBooking) = 0;

    // Common implementation
    //

    G4TNtupleDescription<NT, FT>*  GetNtupleDescriptionInFunction(G4int id,
                                        std::string_view function,
                                        G4bool warn = true) const;
    NT*  GetNtupleInFunction(G4int id,
                             std::string_view function,
                             G4bool warn = true) const;

    // template functions for filling ntuple columns
    template <typename T>
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);

    // Static data members
    static constexpr std::string_view fkClass { "G4TNtupleManager<NT,FT>" };
};

#include "G4TNtupleManager.icc"

#endif

