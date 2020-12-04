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

// The abstract base class for ntuple file output.
//
// Author: Ivana Hrivnacova, 15/09/2020  (ivana@ipno.in2p3.fr)

#ifndef G4VNtupleFileManager_h
#define G4VNtupleFileManager_h 1

#include "G4AnalysisManagerState.hh"
#include "globals.hh"

class G4VNtupleManager;
class G4NtupleBookingManager;

class G4VNtupleFileManager
{
  // Disable using the object managers outside G4VAnalysisManager and
  // its messenger
  friend class G4VAnalysisManager;

  public:
    explicit G4VNtupleFileManager(const G4AnalysisManagerState& state,
                                  const G4String& fileType);
    virtual ~G4VNtupleFileManager();

    // deleted copy constructor & assignment operator
    G4VNtupleFileManager(const G4VNtupleFileManager& rhs) = delete;
    G4VNtupleFileManager& operator=(const G4VNtupleFileManager& rhs) = delete;

    // MT/MPI
    virtual void SetNtupleMerging(G4bool mergeNtuples, 
                   G4int nofReducedNtupleFiles = 0);
    virtual void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true);
    virtual void SetBasketSize(unsigned int basketSize);
    virtual void SetBasketEntries(unsigned int basketEntries);

    virtual std::shared_ptr<G4VNtupleManager> CreateNtupleManager() = 0;
    virtual void SetBookingManager(std::shared_ptr<G4NtupleBookingManager> bookingManager);

    // Methods to be performed at file management
    virtual G4bool ActionAtOpenFile(const G4String& /*fileName*/) = 0;
    virtual G4bool ActionAtWrite() = 0;
    virtual G4bool ActionAtCloseFile(G4bool /*reset*/) = 0;
    virtual G4bool Reset() = 0;
    virtual G4bool IsNtupleMergingSupported() const;

    // Get methods
    G4String GetFileType() const;

  protected:
    // static data members
    const G4AnalysisManagerState& fState;
    G4String  fFileType; 
    std::shared_ptr<G4NtupleBookingManager> fBookingManager;
};

// inline functions

inline void G4VNtupleFileManager::SetBookingManager(
  std::shared_ptr<G4NtupleBookingManager> bookingManager)
{ fBookingManager = bookingManager; }

inline G4bool G4VNtupleFileManager::IsNtupleMergingSupported() const {
  return false;
}

inline G4String G4VNtupleFileManager::GetFileType() const {
  return fFileType;
}

#endif
