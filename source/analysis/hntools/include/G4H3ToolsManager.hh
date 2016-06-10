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
// $Id: G4H3ToolsManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for tools::histo::h3d.
// It implements functions specific to the H3 type
// (defined in g4tools). 
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4H3ToolsManager_h
#define G4H3ToolsManager_h 1

#include "G4VH3Manager.hh"
#include "G4THnManager.hh"
#include "G4HnManager.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>

namespace tools {
namespace histo { 
class h3d; 
}
}

class G4H3ToolsManager : public G4VH3Manager,
                         public G4THnManager<tools::histo::h3d>
{
  public:
    explicit G4H3ToolsManager(const G4AnalysisManagerState& state);
    virtual ~G4H3ToolsManager();

    // Method to add histograms read from a file
    G4int AddH3(const G4String& name, tools::histo::h3d* h3d);
    // Method for merge (MT)
    void AddH3Vector(const std::vector<tools::histo::h3d*>& h3Vector);
    
    // Access methods
    //
    tools::histo::h3d*  GetH3(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    // Iterators
    std::vector<tools::histo::h3d*>::iterator BeginH3();
    std::vector<tools::histo::h3d*>::iterator EndH3();
    std::vector<tools::histo::h3d*>::const_iterator BeginConstH3() const;
    std::vector<tools::histo::h3d*>::const_iterator EndConstH3() const;
                              
    // Access to histogram vector (needed for Write())
    const std::vector<tools::histo::h3d*>& GetH3Vector() const;
    const std::vector<G4HnInformation*>&   GetHnVector() const;
   
  protected:
    // Virtual functions from base class
    //

    // Methods to create histograms
    //
    virtual G4int CreateH3(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4int nzbins, G4double zmin, G4double zmax,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear",
                           const G4String& zbinScheme = "linear") final;
                           
    virtual G4int CreateH3(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") final;
                          
    virtual G4bool SetH3(G4int id,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4int nzbins, G4double zmin, G4double zmax,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear",
                           const G4String& zbinScheme = "linear") final;
                           
    virtual G4bool SetH3(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") final;

    virtual G4bool ScaleH3(G4int id, G4double factor) final;
    
    // Method to fill histograms
    //
    virtual G4bool FillH3(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0) final;
                          

    // Methods to manipulate histograms
    //

    // Access methods
    virtual G4int  GetH3Id(const G4String& name, G4bool warn = true) const final;

    // Access to H3 parameters
    virtual G4int    GetH3Nxbins(G4int id) const final;
    virtual G4double GetH3Xmin(G4int id) const final;
    virtual G4double GetH3Xmax(G4int id) const final;
    virtual G4double GetH3XWidth(G4int id) const final;
    virtual G4int    GetH3Nybins(G4int id) const final;
    virtual G4double GetH3Ymin(G4int id) const final;
    virtual G4double GetH3Ymax(G4int id) const final;
    virtual G4double GetH3YWidth(G4int id) const final;
    virtual G4int    GetH3Nzbins(G4int id) const final;
    virtual G4double GetH3Zmin(G4int id) const final;
    virtual G4double GetH3Zmax(G4int id) const final;
    virtual G4double GetH3ZWidth(G4int id) const final;
        
    // Setters for attributes for plotting
    virtual G4bool SetH3Title(G4int id, const G4String& title) final;
    virtual G4bool SetH3XAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetH3YAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetH3ZAxisTitle(G4int id, const G4String& title) final;

    // Access attributes for plotting
    virtual G4String GetH3Title(G4int id) const final;
    virtual G4String GetH3XAxisTitle(G4int id) const final;
    virtual G4String GetH3YAxisTitle(G4int id) const final;
    virtual G4String GetH3ZAxisTitle(G4int id) const final;
 
     // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output) final;
   
    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() final;

  private:
    void AddH3Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme,
                          G4BinScheme zbinScheme) const;
                            
    // data members
    // static constexpr G4int kDimension = 3;  // not yet supported on vc12
    static const G4int kDimension;
};
// inline methods

inline  std::vector<tools::histo::h3d*>::iterator G4H3ToolsManager::BeginH3()
{ return BeginT(); }

inline  std::vector<tools::histo::h3d*>::iterator G4H3ToolsManager::EndH3()
{ return EndT(); }

inline  std::vector<tools::histo::h3d*>::const_iterator 
G4H3ToolsManager::BeginConstH3() const
{ return BeginConstT(); }

inline  std::vector<tools::histo::h3d*>::const_iterator 
G4H3ToolsManager::EndConstH3() const
{ return EndConstT(); }

inline const std::vector<tools::histo::h3d*>& G4H3ToolsManager::GetH3Vector() const
{ return fTVector; }

inline const std::vector<G4HnInformation*>& G4H3ToolsManager::GetHnVector() const
{ return fHnManager->GetHnVector(); }

inline std::shared_ptr<G4HnManager> G4H3ToolsManager::GetHnManager() 
{ return std::shared_ptr<G4HnManager>(fHnManager); }

#endif

