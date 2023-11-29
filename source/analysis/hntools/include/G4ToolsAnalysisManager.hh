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

// The common implementation of analysis manager classes.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4ToolsAnalysisManager_h
#define G4ToolsAnalysisManager_h 1

#include "G4VAnalysisManager.hh"
#include "G4TH1ToolsManager.hh"
#include "G4TH2ToolsManager.hh"
#include "G4TH3ToolsManager.hh"
#include "G4TP1ToolsManager.hh"
#include "G4TP2ToolsManager.hh"
#include "globals.hh"

#include "G4THnToolsManager.hh"  // make forward declaration if possible

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"

#include <string_view>

class G4PlotManager;

namespace tools {
namespace histo {
class hmpi;
}
}

class G4ToolsAnalysisManager : public  G4VAnalysisManager
{
  friend class G4ToolsAnalysisMessenger;

  public:
    ~G4ToolsAnalysisManager() override;

    // Static methods
    static G4ToolsAnalysisManager* Instance();
    static G4bool IsInstance();

    // Access methods
    tools::histo::h1d*  GetH1(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::h2d*  GetH2(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::h3d*  GetH3(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::p1d*  GetP1(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::p2d*  GetP2(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;

    // Iterators
    std::vector<tools::histo::h1d*>::iterator BeginH1();
    std::vector<tools::histo::h1d*>::iterator EndH1();
    std::vector<tools::histo::h1d*>::const_iterator BeginConstH1() const;
    std::vector<tools::histo::h1d*>::const_iterator EndConstH1() const;

    std::vector<tools::histo::h2d*>::iterator BeginH2();
    std::vector<tools::histo::h2d*>::iterator EndH2();
    std::vector<tools::histo::h2d*>::const_iterator BeginConstH2() const;
    std::vector<tools::histo::h2d*>::const_iterator EndConstH2() const;

    std::vector<tools::histo::h3d*>::iterator BeginH3();
    std::vector<tools::histo::h3d*>::iterator EndH3();
    std::vector<tools::histo::h3d*>::const_iterator BeginConstH3() const;
    std::vector<tools::histo::h3d*>::const_iterator EndConstH3() const;

    std::vector<tools::histo::p1d*>::iterator BeginP1();
    std::vector<tools::histo::p1d*>::iterator EndP1();
    std::vector<tools::histo::p1d*>::const_iterator BeginConstP1() const;
    std::vector<tools::histo::p1d*>::const_iterator EndConstP1() const;

    std::vector<tools::histo::p2d*>::iterator BeginP2();
    std::vector<tools::histo::p2d*>::iterator EndP2();
    std::vector<tools::histo::p2d*>::const_iterator BeginConstP2() const;
    std::vector<tools::histo::p2d*>::const_iterator EndConstP2() const;

  protected:
    explicit G4ToolsAnalysisManager(const G4String& type);

    // Virtual methods from base class
    G4bool OpenFileImpl(const G4String& fileName) override;
    G4bool WriteImpl() override;
    G4bool CloseFileImpl(G4bool reset) override;
    G4bool ResetImpl() override;
    void ClearImpl() override;
    G4bool PlotImpl() final;
    G4bool MergeImpl(tools::histo::hmpi* hmpi) final;
    G4bool IsOpenFileImpl() const final;

    // Methods
    G4bool IsEmpty();

     // Static data members
    inline static G4ToolsAnalysisManager* fgMasterToolsInstance { nullptr };
    inline static G4ThreadLocal G4ToolsAnalysisManager* fgToolsInstance { nullptr };
    static constexpr std::string_view fkClass { "G4ToolsAnalysisManager" };

    // Data members
    G4THnToolsManager<kDim1, tools::histo::h1d>* fH1Manager { nullptr };
    G4THnToolsManager<kDim2, tools::histo::h2d>* fH2Manager { nullptr };
    G4THnToolsManager<kDim3, tools::histo::h3d>* fH3Manager { nullptr };
    G4THnToolsManager<kDim2, tools::histo::p1d>* fP1Manager { nullptr };
    G4THnToolsManager<kDim3, tools::histo::p2d>* fP2Manager { nullptr };

  private:
    //  // Static data members
    // Static G4ThreadLocal G4ToolsAnalysisManager* fgToolsInstance;
    // Methods
    template <typename HT>
    G4bool WriteT(const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector);

    G4bool WriteHns();
    G4bool ResetHns();
    G4bool MergeHns();

    // Data members
    std::shared_ptr<G4PlotManager>   fPlotManager { nullptr };
 };

#include "G4ToolsAnalysisManager.icc"

#endif
