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

// Base class for readers file managers.
//
// Author: Ivana Hrivnacova, 06/08/2021  (ivana@ijclab.in2p3.fr)

#ifndef G4VRFileManager_h
#define G4VRFileManager_h 1

#include "G4BaseFileManager.hh"
#include "G4VTHnRFileManager.hh"

#include <memory>
#include <string_view>

namespace tools {
namespace histo {
class h1d;
class h2d;
class h3d;
class p1d;
class p2d;
}
}

class G4VRFileManager : public G4BaseFileManager
{
  public:
    explicit G4VRFileManager(const G4AnalysisManagerState& state)
      : G4BaseFileManager(state) {}
    G4VRFileManager() = delete;
    ~G4VRFileManager() override = default;

    // Methods applied to all registered files
    virtual void CloseFiles() = 0;

    // Access to helpers
    template <typename HT>
    std::shared_ptr<G4VTHnRFileManager<HT>> GetHnRFileManager() const;

  protected:
    // Static data members
    static constexpr std::string_view fkClass { "G4VRFileManager" };

    // Data members
    // FileManagers per object type
    std::shared_ptr<G4VTHnRFileManager<tools::histo::h1d>> fH1RFileManager { nullptr };
    std::shared_ptr<G4VTHnRFileManager<tools::histo::h2d>> fH2RFileManager { nullptr };
    std::shared_ptr<G4VTHnRFileManager<tools::histo::h3d>> fH3RFileManager { nullptr };
    std::shared_ptr<G4VTHnRFileManager<tools::histo::p1d>> fP1RFileManager { nullptr };
    std::shared_ptr<G4VTHnRFileManager<tools::histo::p2d>> fP2RFileManager { nullptr };
};

// inline functions

template <>
inline
std::shared_ptr<G4VTHnRFileManager<tools::histo::h1d>>
G4VRFileManager::GetHnRFileManager<tools::histo::h1d>() const
{ return fH1RFileManager; }

template <>
inline
std::shared_ptr<G4VTHnRFileManager<tools::histo::h2d>>
G4VRFileManager::GetHnRFileManager<tools::histo::h2d>() const
{ return fH2RFileManager; }

template <>
inline
std::shared_ptr<G4VTHnRFileManager<tools::histo::h3d>>
G4VRFileManager::GetHnRFileManager<tools::histo::h3d>() const
{ return fH3RFileManager; }

template <>
inline
std::shared_ptr<G4VTHnRFileManager<tools::histo::p1d>>
G4VRFileManager::GetHnRFileManager<tools::histo::p1d>() const
{ return fP1RFileManager; }

template <>
inline
std::shared_ptr<G4VTHnRFileManager<tools::histo::p2d>>
G4VRFileManager::GetHnRFileManager<tools::histo::p2d>() const
{ return fP2RFileManager; }

#endif
