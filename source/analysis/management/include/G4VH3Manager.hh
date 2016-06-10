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
// $Id: G4VH3Manager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Base class for H3 manager.
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4VH3Manager_h
#define G4VH3Manager_h 1

#include "globals.hh"

#include <vector>
#include <memory>

class G4HnManager;

class G4VH3Manager
{
  // Disable using the object managers outside G4VAnalysisManager
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;

  public:
    G4VH3Manager() {}
    virtual ~G4VH3Manager() {}

    // deleted copy constructor & assignment operator
    G4VH3Manager(const G4VH3Manager& rhs) = delete;
    G4VH3Manager& operator=(const G4VH3Manager& rhs) = delete;

  protected:
    // Methods for handling histograms
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
                           const G4String& zbinScheme = "linear") = 0;
                           
    virtual G4int CreateH3(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") = 0;
                           
    virtual G4bool SetH3(G4int id,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nzbins, G4double zmin, G4double zmax,
                           G4int nybins, G4double ymin, G4double ymax,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear",
                           const G4String& zbinScheme = "linear") = 0;

    virtual G4bool SetH3(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const std::vector<G4double>& zedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") = 0;

    virtual G4bool ScaleH3(G4int id, G4double factor) = 0;

    // Methods to fill histograms
    virtual G4bool FillH3(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0) = 0;
    
    // Access methods
    virtual G4int  GetH3Id(const G4String& name, G4bool warn = true) const = 0;

    // Access to H3 parameters
    virtual G4int    GetH3Nxbins(G4int id) const = 0;
    virtual G4double GetH3Xmin(G4int id) const = 0;
    virtual G4double GetH3Xmax(G4int id) const = 0;
    virtual G4double GetH3XWidth(G4int id) const = 0;
    virtual G4int    GetH3Nybins(G4int id) const = 0;
    virtual G4double GetH3Ymin(G4int id) const = 0;
    virtual G4double GetH3Ymax(G4int id) const = 0;
    virtual G4double GetH3YWidth(G4int id) const = 0;
    virtual G4int    GetH3Nzbins(G4int id) const = 0;
    virtual G4double GetH3Zmin(G4int id) const = 0;
    virtual G4double GetH3Zmax(G4int id) const = 0;
    virtual G4double GetH3ZWidth(G4int id) const = 0;

    // Setters for attributes for plotting
    virtual G4bool SetH3Title(G4int id, const G4String& title) = 0;
    virtual G4bool SetH3XAxisTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetH3YAxisTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetH3ZAxisTitle(G4int id, const G4String& title) = 0;

    // Access attributes for plotting
    virtual G4String GetH3Title(G4int id) const = 0;
    virtual G4String GetH3XAxisTitle(G4int id) const = 0;
    virtual G4String GetH3YAxisTitle(G4int id) const = 0;
    virtual G4String GetH3ZAxisTitle(G4int id) const = 0;

    // Methods to manipulate histograms
    virtual G4bool WriteOnAscii(std::ofstream& output) = 0;

    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() = 0;
};

#endif

