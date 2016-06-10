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
// $Id: G4H1ToolsManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for tools::histo::h1d.
// It implements functions specific to the H1 type
// (defined in g4tools). 
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4H1ToolsManager_h
#define G4H1ToolsManager_h 1

#include "G4VH1Manager.hh"
#include "G4THnManager.hh"
#include "G4HnManager.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>

namespace tools {
namespace histo { 
class h1d; 
}
}

class G4H1ToolsManager : public G4VH1Manager,
                         public G4THnManager<tools::histo::h1d>
{
  public:
    explicit G4H1ToolsManager(const G4AnalysisManagerState& state);
    virtual ~G4H1ToolsManager();

    // Method to add histograms read from a file
    G4int AddH1(const G4String& name, tools::histo::h1d* h1d);
    // Method for merge (MT)
    void  AddH1Vector(const std::vector<tools::histo::h1d*>& h1Vector);
    
    // Access methods
    //
    tools::histo::h1d*  GetH1(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
                              
    // Iterators
    std::vector<tools::histo::h1d*>::iterator BeginH1();
    std::vector<tools::histo::h1d*>::iterator EndH1();
    std::vector<tools::histo::h1d*>::const_iterator BeginConstH1() const;
    std::vector<tools::histo::h1d*>::const_iterator EndConstH1() const;
                              
    // Access to histogram vector (needed for Write())
    const std::vector<tools::histo::h1d*>& GetH1Vector() const;
    const std::vector<G4HnInformation*>&   GetHnVector() const;  
    
  protected:
    // Virtual functions from base class
    //

    // Methods to create histograms
    //
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binScheme = "linear") final;
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none") final;
                           
    virtual G4bool SetH1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binSchemeName = "linear") final;
    virtual G4bool SetH1(G4int id,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none") final;
    virtual G4bool ScaleH1(G4int id, G4double factor) final;
    
    // Method to fill histograms
    //
    virtual G4bool FillH1(G4int id, G4double value, G4double weight = 1.0) final;

    // Access methods
    //
    virtual G4int  GetH1Id(const G4String& name, G4bool warn = true) const final;

    // Access to H1 parameters
    virtual G4int    GetH1Nbins(G4int id) const final;
    virtual G4double GetH1Xmin(G4int id) const final;
    virtual G4double GetH1Xmax(G4int id) const final;
    virtual G4double GetH1Width(G4int id) const final;

    // Attributes for plotting
    //

    // Setters
    virtual G4bool SetH1Title(G4int id, const G4String& title) final;
    virtual G4bool SetH1XAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetH1YAxisTitle(G4int id, const G4String& title) final;

    // Accessors
    virtual G4String GetH1Title(G4int id) const final;
    virtual G4String GetH1XAxisTitle(G4int id) const final;
    virtual G4String GetH1YAxisTitle(G4int id) const final;

    // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output) final;

    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() final;

  private:
    // methods
    //
    void AddH1Information(const G4String& name,  
                          const G4String& unitName, 
                          const G4String& fcnName,
                          G4BinScheme binScheme) const;

    // data members
    //static constexpr G4int kDimension = 1;  // not yet supported on vc12
    static const G4int kDimension;
};

// inline methods

inline  std::vector<tools::histo::h1d*>::iterator G4H1ToolsManager::BeginH1()
{ return BeginT(); }

inline  std::vector<tools::histo::h1d*>::iterator G4H1ToolsManager::EndH1()
{ return EndT(); }

inline  std::vector<tools::histo::h1d*>::const_iterator 
G4H1ToolsManager::BeginConstH1() const
{ return BeginConstT(); }

inline  std::vector<tools::histo::h1d*>::const_iterator 
G4H1ToolsManager::EndConstH1() const
{ return EndConstT(); }

inline const std::vector<tools::histo::h1d*>& G4H1ToolsManager::GetH1Vector() const
{ return fTVector; }

inline const std::vector<G4HnInformation*>& G4H1ToolsManager::GetHnVector() const
{ return fHnManager->GetHnVector(); }

inline std::shared_ptr<G4HnManager> G4H1ToolsManager::GetHnManager() 
{ return std::shared_ptr<G4HnManager>(fHnManager); }

#endif

