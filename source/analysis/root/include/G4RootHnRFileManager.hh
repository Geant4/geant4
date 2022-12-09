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

// The manager class for histogram/profile reading from a Root file.

// Author: Ivana Hrivnacova, 26/08/2021  (ivana@ipno.in2p3.fr)

#ifndef G4RootHnRFileManager_h
#define G4RootHnRFileManager_h 1

#include "G4VTHnRFileManager.hh"

#include <string_view>
#include <tuple>

class G4RootRFileManager;

namespace tools {
namespace rroot {
class buffer;
class TDirectory;
}
}

template <typename HT>
class G4RootHnRFileManager : public G4VTHnRFileManager<HT>
{
  public:
    explicit G4RootHnRFileManager(G4RootRFileManager* rfileManger)
      : G4VTHnRFileManager<HT>(), fRFileManager(rfileManger) {}
    G4RootHnRFileManager() = delete;
    ~G4RootHnRFileManager() override = default;

    // Methods for reading objects
    HT* Read(const G4String& htName, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;

  private:
    // Methods
    tools::rroot::buffer* GetBuffer(
                             const G4String& fileName,
                             const G4String& dirName,
                             const G4String& name);
    HT* ReadT(tools::rroot::buffer* buffer);

    // Static data members
    static constexpr std::string_view fkClass { "G4RootHnRFileManager<HT>" };
    // Data members
    G4RootRFileManager* fRFileManager { nullptr };
};

// inline functions

#include "G4RootHnRFileManager.icc"

#endif
