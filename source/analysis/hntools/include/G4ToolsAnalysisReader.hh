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

// The common implementation of analysis reader classes.

// Author: Ivana Hrivnacova, 20/07/2015  (ivana@ipno.in2p3.fr)

#ifndef G4ToolsAnalysisReader_h
#define G4ToolsAnalysisReader_h 1

#include "G4VAnalysisReader.hh"
#include "G4THnManager.hh"
#include "globals.hh"

#include "G4THnToolsManager.hh"  // TO DO: make forward declaration if possible

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"

class G4P2ToolsManager;

namespace tools {
namespace histo {
class hmpi;
}
}

class G4ToolsAnalysisReader : public  G4VAnalysisReader
{
  public:
    ~G4ToolsAnalysisReader() override = default;

    // Access methods
    tools::histo::h1d*  GetH1(G4int id, G4bool warn = true) const;
    tools::histo::h2d*  GetH2(G4int id, G4bool warn = true) const;
    tools::histo::h3d*  GetH3(G4int id, G4bool warn = true) const;
    tools::histo::p1d*  GetP1(G4int id, G4bool warn = true) const;
    tools::histo::p2d*  GetP2(G4int id, G4bool warn = true) const;

  protected:
    explicit G4ToolsAnalysisReader(const G4String& type);

    // Virtual methods from base class
    G4int ReadH1Impl(const G4String& h1Name, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;
    G4int ReadH2Impl(const G4String& h2Name, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;
    G4int ReadH3Impl(const G4String& h3Name, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;
    G4int ReadP1Impl(const G4String& p1Name, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;
    G4int ReadP2Impl(const G4String& p2Name, const G4String& fileName, const G4String& dirName,
      G4bool isUserFileName) final;

    // Fuction specific to output type
    template <typename HT>
    G4int ReadTImpl(const G4String& htName, const G4String& fileName,
                    const G4String& dirName, G4bool isUserFileName,
                    G4THnManager<HT>* htManager);

    // Methods
    G4bool Reset();

    // Static data members
    static constexpr std::string_view fkClass { "G4ToolsAnalysisReader" };

    // Data members
    G4THnToolsManager<1, tools::histo::h1d>* fH1Manager { nullptr };
    G4THnToolsManager<2, tools::histo::h2d>* fH2Manager { nullptr };
    G4THnToolsManager<3, tools::histo::h3d>* fH3Manager { nullptr };
    G4THnToolsManager<2, tools::histo::p1d>* fP1Manager { nullptr };
    G4THnToolsManager<3, tools::histo::p2d>* fP2Manager { nullptr };
 };

#include "G4ToolsAnalysisReader.icc"

#endif
