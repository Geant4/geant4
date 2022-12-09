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

// Manager class for ntuple booking for all output types.
//
// Author: Ivana Hrivnacova, 19/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4NtupleBookingManager_h
#define G4NtupleBookingManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include "tools/ntuple_booking"

#include <vector>
#include <string_view>

struct G4NtupleBooking
{
  G4NtupleBooking() = default;
  ~G4NtupleBooking() = default;

  tools::ntuple_booking fNtupleBooking;
  G4int fNtupleId { G4Analysis::kInvalidId };
  G4String fFileName;
  G4bool fActivation { true };
};

class G4NtupleBookingManager : public G4BaseAnalysisManager
{
  // Disable using the object managers outside G4VAnalysisManager and
  // its messenger
  friend class G4VAnalysisManager;

  public:
    explicit G4NtupleBookingManager(const G4AnalysisManagerState& state);
    G4NtupleBookingManager() = delete;
    ~G4NtupleBookingManager() override;

    const std::vector<G4NtupleBooking*>& GetNtupleBookingVector() const;

    // File type
    // (only one output file type allowed for ntuples)
    void SetFileType(const G4String& fileType);
    G4String GetFileType() const;
    G4bool IsEmpty() const;

  protected:
    // Methods to create ntuples
    //
    G4int CreateNtuple(const G4String& name, const G4String& title);

    // Create columns in the last created ntuple (from base class)
    // Create columns in the last created ntuple
    G4int CreateNtupleIColumn(const G4String& name,
                              std::vector<int>* vector);
    G4int CreateNtupleFColumn(const G4String& name,
                              std::vector<float>* vector) ;
    G4int CreateNtupleDColumn(const G4String& name,
                              std::vector<double>* vector);
    G4int CreateNtupleSColumn(const G4String& name,
                              std::vector<std::string>* vector);
    G4NtupleBooking*  FinishNtuple() ;

    // Create columns in the ntuple with given id
    G4int CreateNtupleIColumn(G4int ntupleId,
                    const G4String& name, std::vector<int>* vector);
    G4int CreateNtupleFColumn(G4int ntupleId,
                    const G4String& name, std::vector<float>* vector);
    G4int CreateNtupleDColumn(G4int ntupleId,
                    const G4String& name, std::vector<double>* vector);
    G4int CreateNtupleSColumn(G4int ntupleId,
                    const G4String& name, std::vector<std::string>* vector);
    G4NtupleBooking*  FinishNtuple(G4int ntupleId);

    // The ntuple column ids are generated automatically starting from 0;
    // with the following function it is possible to change it
    // to start from another value
    G4bool SetFirstNtupleColumnId(G4int firstId);
    G4int  GetFirstNtupleColumnId() const;

    // Activation option
    //
    void  SetActivation(G4bool activation);
    void  SetActivation(G4int ntupleId, G4bool activation);
    G4bool  GetActivation(G4int ntupleId) const;

    // // File name option
    void  SetFileName(const G4String& fileName);
    void  SetFileName(G4int id, const G4String& fileName);
    G4String GetFileName(G4int id) const;

    // Access methods
    G4int GetNofNtupleBookings() const;
    // G4int GetNofNtuples() const;

    // Clear data
    void ClearData();

    // Data members
    std::vector<G4NtupleBooking*> fNtupleBookingVector;

  private:
    // Methods
    // Common implementation
    G4NtupleBooking*  GetNtupleBookingInFunction(
                                        G4int id,
                                        std::string_view function,
                                        G4bool warn = true) const;

    G4bool CheckName(const G4String& name, const G4String& objectType) const;
    template <typename T>
    G4int CreateNtupleTColumn(G4int ntupleId,
                    const G4String& name, std::vector<T>* vector);
    G4int GetCurrentNtupleId() const;

    // Static data members
    static constexpr std::string_view fkClass { "G4NtupleBookingManager" };

    // Data members
    G4String fFileType;
    G4int    fFirstNtupleColumnId { 0 };
    G4bool   fLockFirstNtupleColumnId { false };
};

#include "G4NtupleBookingManager.icc"

#endif
