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

// The manager for histogram/profile Hfd5 file output.

// Author: Ivana Hrivnacova, 15/09/2020  (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5HnFileManager_h
#define G4Hdf5HnFileManager_h 1

#include "G4VTHnFileManager.hh"

#include "toolx/hdf5/ntuple" // for hid_t

#include <string_view>

class G4Hdf5FileManager;

template <typename HT>
class G4Hdf5HnFileManager : public G4VTHnFileManager<HT>
{
  public:
    explicit G4Hdf5HnFileManager(G4Hdf5FileManager* fileManger)
      : G4VTHnFileManager<HT>(), fFileManager(fileManger) {}
    G4Hdf5HnFileManager() = delete;
    ~G4Hdf5HnFileManager() override = default;

    // Methods for writing objects
    // Write to a new file (the file is closed after write)
    G4bool WriteExtra(HT* ht, const G4String& htName, const G4String& fileName) final;
    // Write to the default file  (handled with OpenFile()/CloseFile methods)
    G4bool Write(HT* ht, const G4String& htName, G4String& fileName) final;

  private:
    // Methods
    G4bool WriteImpl(hid_t hdirectory, HT* ht, const G4String& htName);
    // Static data members
    static constexpr std::string_view fkClass { "G4Hdf5HnFileManager<HT>" };
    // Data members
    G4Hdf5FileManager* fFileManager { nullptr };
};

// inline functions

#include "G4Hdf5HnFileManager.icc"

#endif
