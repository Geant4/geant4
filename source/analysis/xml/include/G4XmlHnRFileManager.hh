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

// The manager class for histogram/profile reading from a Xml file.

// Author: Ivana Hrivnacova, 26/08/2021  (ivana@ipno.in2p3.fr)

#ifndef G4XmlHnRFileManager_h
#define G4XmlHnRFileManager_h 1

#include "G4VTHnRFileManager.hh"

#include <string_view>

class G4XmlRFileManager;

template <typename HT>
class G4XmlHnRFileManager : public G4VTHnRFileManager<HT>
{
  public:
    G4XmlHnRFileManager(G4XmlRFileManager* rfileManger)
      : G4VTHnRFileManager<HT>(), fRFileManager(rfileManger) {}
    G4XmlHnRFileManager() = delete;
    ~G4XmlHnRFileManager() override = default;

    // Methods for writing objects
    HT* Read(const G4String& htName, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4XmlHnRFileManager<HT>" };
    // Data members
    G4XmlRFileManager* fRFileManager { nullptr };
};

// inline functions

#include "G4XmlHnRFileManager.icc"

#endif
