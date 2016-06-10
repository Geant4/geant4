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

// Manager class for tools::histo::P1d.
// It implements functions specific to the P1 type
// (defined in g4tools). 
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4P1ToolsManager_h
#define G4P1ToolsManager_h 1

#include "G4VP1Manager.hh"
#include "G4THnManager.hh"
#include "G4HnManager.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>

namespace tools {
namespace histo { 
class p1d; 
}
}

class G4P1ToolsManager : public G4VP1Manager,
                         public G4THnManager<tools::histo::p1d>
{
  public:
    explicit G4P1ToolsManager(const G4AnalysisManagerState& state);
    virtual ~G4P1ToolsManager();

    // Method to add profiles read from a file
    G4int AddP1(const G4String& name, tools::histo::p1d* p1d);
    // Method for merge (MT)
    void  AddP1Vector(const std::vector<tools::histo::p1d*>& p1Vector);
    
    // Access methods
    //
    tools::histo::p1d*  GetP1(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
                              
    // Iterators
    std::vector<tools::histo::p1d*>::iterator BeginP1();
    std::vector<tools::histo::p1d*>::iterator EndP1();
    std::vector<tools::histo::p1d*>::const_iterator BeginConstP1() const;
    std::vector<tools::histo::p1d*>::const_iterator EndConstP1() const;
                              
    // Access to profile vector (needed for Write())
    const std::vector<tools::histo::p1d*>& GetP1Vector() const;
    const std::vector<G4HnInformation*>&   GetHnVector() const;  
    
  protected:
    // Virtual functions from base class
    //

    // Methods to create profiles
    //
    virtual G4int CreateP1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& xbinScheme = "linear") final;
    virtual G4int CreateP1(const G4String& name, const G4String& title,
                           const std::vector<G4double>& edges,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none") final;
                           
    virtual G4bool SetP1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& xbinScheme = "linear") final;
    virtual G4bool SetP1(G4int id,
                           const std::vector<G4double>& edges,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none") final;
    virtual G4bool ScaleP1(G4int id, G4double factor) final;
    
    // Method to fill profiles
    //
    virtual G4bool FillP1(G4int id, G4double xvalue, G4double yvalue, 
                          G4double weight = 1.0) final;

    // Access methods
    //
    virtual G4int  GetP1Id(const G4String& name, G4bool warn = true) const final;

    // Access to P1 parameters
    virtual G4int    GetP1Nbins(G4int id) const final;
    virtual G4double GetP1Xmin(G4int id) const final;
    virtual G4double GetP1Xmax(G4int id) const final;
    virtual G4double GetP1XWidth(G4int id) const final;
    virtual G4double GetP1Ymin(G4int id) const final;
    virtual G4double GetP1Ymax(G4int id) const final;

    // Attributes for plotting
    //

    // Setters
    virtual G4bool SetP1Title(G4int id, const G4String& title) final;
    virtual G4bool SetP1XAxisTitle(G4int id, const G4String& title) final;
    virtual G4bool SetP1YAxisTitle(G4int id, const G4String& title) final;

    // Accessors
    virtual G4String GetP1Title(G4int id) const final;
    virtual G4String GetP1XAxisTitle(G4int id) const final;
    virtual G4String GetP1YAxisTitle(G4int id) const final;

    // Write data on ASCII file
    //virtual G4bool WriteOnAscii(std::ofstream& output) final;

    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() final;

  private:
    // methods
    //
    void AddP1Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme) const;

    // data members
    //static constexpr G4int kDimension = 1;  // not yet supported on vc12
    static const G4int kDimension;
};

// inline methods

inline  std::vector<tools::histo::p1d*>::iterator G4P1ToolsManager::BeginP1()
{ return BeginT(); }

inline  std::vector<tools::histo::p1d*>::iterator G4P1ToolsManager::EndP1()
{ return EndT(); }

inline  std::vector<tools::histo::p1d*>::const_iterator 
G4P1ToolsManager::BeginConstP1() const
{ return BeginConstT(); }

inline  std::vector<tools::histo::p1d*>::const_iterator 
G4P1ToolsManager::EndConstP1() const
{ return EndConstT(); }

inline const std::vector<tools::histo::p1d*>& G4P1ToolsManager::GetP1Vector() const
{ return fTVector; }

inline const std::vector<G4HnInformation*>& G4P1ToolsManager::GetHnVector() const
{ return fHnManager->GetHnVector(); }

inline std::shared_ptr<G4HnManager> G4P1ToolsManager::GetHnManager()
{ return std::shared_ptr<G4HnManager>(fHnManager); }

#endif

