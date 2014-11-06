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
// $Id: ExG4HbookH3DummyManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

/// \file hbook/include/ExG4HbookH3DummyManager.hh
/// \brief Definition of the ExG4HbookH3DummyManager class

// Author: Ivana Hrivnacova, 03/11/2014  (ivana@ipno.in2p3.fr)

#ifndef ExG4HbookH3DummyManager_h
#define ExG4HbookH3DummyManager_h 1

#include "G4VH3Manager.hh"
#include "G4BaseToolsManager.hh"
#include "G4HnManager.hh"
#include "globals.hh"

#include <vector>
#include <map>

namespace tools {
namespace histo { 
class h3d; 
}
}

/// Manager class for H3 with dummy implementation. 
///
/// It will just issue warnings. 

class ExG4HbookH3DummyManager : public G4VH3Manager
{
  public:
    ExG4HbookH3DummyManager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookH3DummyManager();
   
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
                           const G4String& zbinScheme = "linear");
                           
    virtual G4int CreateH3(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none");
                          
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
                           const G4String& zbinScheme = "linear");
                           
    virtual G4bool SetH3(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none");

    virtual G4bool ScaleH3(G4int id, G4double factor);
    
    // Method to fill histograms
    //
    virtual G4bool FillH3(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0);
                          

    // Methods to manipulate histograms
    //

    // Access methods
    virtual G4int  GetH3Id(const G4String& name, G4bool warn = true) const;

    // Access to H3 parameters
    virtual G4int    GetH3Nxbins(G4int id) const;
    virtual G4double GetH3Xmin(G4int id) const;
    virtual G4double GetH3Xmax(G4int id) const;
    virtual G4double GetH3XWidth(G4int id) const;
    virtual G4int    GetH3Nybins(G4int id) const;
    virtual G4double GetH3Ymin(G4int id) const;
    virtual G4double GetH3Ymax(G4int id) const;
    virtual G4double GetH3YWidth(G4int id) const;
    virtual G4int    GetH3Nzbins(G4int id) const;
    virtual G4double GetH3Zmin(G4int id) const;
    virtual G4double GetH3Zmax(G4int id) const;
    virtual G4double GetH3ZWidth(G4int id) const;
        
    // Setters for attributes for plotting
    virtual G4bool SetH3Title(G4int id, const G4String& title);
    virtual G4bool SetH3XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH3YAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH3ZAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    virtual G4String GetH3Title(G4int id) const;
    virtual G4String GetH3XAxisTitle(G4int id) const;
    virtual G4String GetH3YAxisTitle(G4int id) const;
    virtual G4String GetH3ZAxisTitle(G4int id) const;
 
     // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output);
    
  private:
    // methods
    void ExceptionForHistograms(const G4String& functionName);
    void ExceptionForHistogramsConst(const G4String& functionName) const;

    // data members
    G4bool fWarn;    
};

#endif

