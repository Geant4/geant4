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

// The manager class for MPI applications.

// Author: Ivana Hrivnacova, 25/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4MPIToolsManager_h
#define G4MPIToolsManager_h 1

#include "G4AnalysisManagerState.hh"
#include "G4HnInformation.hh"
#include "G4ios.hh"

#include "tools/impi_world"
#include "tools/histo/hmpi"

#include <vector>
#include <string_view>

class G4MPIToolsManager
{
  public:
    G4MPIToolsManager(const G4AnalysisManagerState& state,
                      tools::histo::hmpi* hmpi)
    : fState(state), fHmpi(hmpi) {}
    G4MPIToolsManager() = delete;
    virtual ~G4MPIToolsManager() = default;

  public:
    // Methods
    template <typename HT>
    G4bool Merge(const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector

);
  private:
    // Methods
    template <typename HT>
    G4bool Send(G4int nofActiveT,
                const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector

);

    template <typename HT>
    G4bool Receive(G4int nofActiveT,
                const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector

);

    // Static data members
    static constexpr std::string_view fkClass { "G4MPIToolsManager" };

    // Data members
    const G4AnalysisManagerState& fState;
    tools::histo::hmpi*  fHmpi;
};

// inline functions
#include "G4MPIToolsManager.icc"

#endif
