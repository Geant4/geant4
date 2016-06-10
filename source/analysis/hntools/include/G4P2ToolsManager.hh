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
// $Id$

// Manager class for tools::histo::p2d.
// It implements functions specific to the P2 type
// (defined in g4tools). 
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4P2ToolsManager_h
#define G4P2ToolsManager_h 1

#include "G4VP2Manager.hh"
#include "G4THnManager.hh"
#include "G4HnManager.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>

namespace tools {
namespace histo { 
class p2d; 
}
}

class G4P2ToolsManager : public G4VP2Manager,
                         public G4THnManager<tools::histo::p2d>
{
  public:
    explicit G4P2ToolsManager(const G4AnalysisManagerState& state);
    virtual ~G4P2ToolsManager();

    // Method to add histograms read from a file
    G4int AddP2(const G4String& name, tools::histo::p2d* p2d);
    // Method for merge (MT)
    void AddP2Vector(const std::vector<tools::histo::p2d*>& p2Vector);
    
    // Access methods
    //
    tools::histo::p2d*  GetP2(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    // Iterators
    std::vector<tools::histo::p2d*>::iterator BeginP2();
    std::vector<tools::histo::p2d*>::iterator EndP2();
    std::vector<tools::histo::p2d*>::const_iterator BeginConstP2() const;
    std::vector<tools::histo::p2d*>::const_iterator EndConstP2() const;
                              
    // Access to histogram vector (needed for Write())
    const std::vector<tools::histo::p2d*>& GetP2Vector() const;
    const std::vector<G4HnInformation*>&   GetHnVector() const;
   
  protected:
    // Virtual functions from base class
    //

    // Methods to create histograms
    //
    virtual G4int CreateP2(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear") final;
                           
    virtual G4int CreateP2(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") final;
                          
    virtual G4bool SetP2(G4int id,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear") final;
                           
    virtual G4bool SetP2(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") final;

    virtual G4bool ScaleP2(G4int id, G4double factor) final;
    
    // Method to fill histograms
    //
    virtual G4bool FillP2(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0) final;
                          

    // Methods to manipulate histograms
    //

    // Access methods
    virtual G4int  GetP2Id(const G4String& name, G4bool warn = true) const final;

    // Access to P2 parameters
    virtual G4int    GetP2Nxbins(G4int id) const final;
    virtual G4double GetP2Xmin(G4int id) const final;
    virtual G4double GetP2Xmax(G4int id) const final;
    virtual G4double GetP2XWidth(G4int id) const final;
    virtual G4int    GetP2Nybins(G4int id) const final;
    virtual G4double GetP2Ymin(G4int id) const final;
    virtual G4double GetP2Ymax(G4int id) const final;
    virtual G4double GetP2YWidth(G4int id) const final;
    virtual G4double GetP2Zmin(G4int id) const final;
    virtual G4double GetP2Zmax(G4int id) const final;
        
    // Setters for attributes for plotting
    virtual G4bool SetP2Title(G4int id, const G4String& title) final;
    virtual G4bool SetP2XAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetP2YAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetP2ZAxisTitle(G4int id, const G4String& title) final;

    // Access attributes for plotting
    virtual G4String GetP2Title(G4int id) const final;
    virtual G4String GetP2XAxisTitle(G4int id) const final;
    virtual G4String GetP2YAxisTitle(G4int id) const final;
    virtual G4String GetP2ZAxisTitle(G4int id) const final;
 
     // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output) final;
   
    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() final;

  private:
    void AddP2Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme) const;
                            
    // data members
    // static constexpr G4int kDimension = 2;  // not yet supported on vc12
    static const G4int kDimension;
};

// inline methods

inline  std::vector<tools::histo::p2d*>::iterator G4P2ToolsManager::BeginP2()
{ return BeginT(); }

inline  std::vector<tools::histo::p2d*>::iterator G4P2ToolsManager::EndP2()
{ return EndT(); }

inline  std::vector<tools::histo::p2d*>::const_iterator 
G4P2ToolsManager::BeginConstP2() const
{ return BeginConstT(); }

inline  std::vector<tools::histo::p2d*>::const_iterator 
G4P2ToolsManager::EndConstP2() const
{ return EndConstT(); }

inline const std::vector<tools::histo::p2d*>& G4P2ToolsManager::GetP2Vector() const
{ return fTVector; }

inline const std::vector<G4HnInformation*>& G4P2ToolsManager::GetHnVector() const
{ return fHnManager->GetHnVector(); }

inline std::shared_ptr<G4HnManager> G4P2ToolsManager::GetHnManager() 
{ return std::shared_ptr<G4HnManager>(fHnManager); }

#endif

