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

// Base class for P2 manager.
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4VP2Manager_h
#define G4VP2Manager_h 1

#include "globals.hh"

#include <vector>
#include <memory>

class G4HnManager;

class G4VP2Manager
{
  // Disable using the object managers outside G4VAnalysisManager
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;

  public:
    G4VP2Manager() {}
    virtual ~G4VP2Manager() {}

    // deleted copy constructor & assignment operator
    G4VP2Manager(const G4VP2Manager& rhs) = delete;
    G4VP2Manager& operator=(const G4VP2Manager& rhs) = delete;

  protected:
    // Methods for handling histograms
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
                           const G4String& ybinScheme = "linear") = 0;
                           
    virtual G4int CreateP2(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") = 0;
                           
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
                           const G4String& ybinScheme = "linear") = 0;

    virtual G4bool SetP2(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none") = 0;

    virtual G4bool ScaleP2(G4int id, G4double factor) = 0;

    // Methods to fill histograms
    virtual G4bool FillP2(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0) = 0;
    
    // Access methods
    virtual G4int  GetP2Id(const G4String& name, G4bool warn = true) const = 0;

    // Access to P2 parameters
    virtual G4int    GetP2Nxbins(G4int id) const = 0;
    virtual G4double GetP2Xmin(G4int id) const = 0;
    virtual G4double GetP2Xmax(G4int id) const = 0;
    virtual G4double GetP2XWidth(G4int id) const = 0;
    virtual G4int    GetP2Nybins(G4int id) const = 0;
    virtual G4double GetP2Ymin(G4int id) const = 0;
    virtual G4double GetP2Ymax(G4int id) const = 0;
    virtual G4double GetP2YWidth(G4int id) const = 0;
    virtual G4double GetP2Zmin(G4int id) const = 0;
    virtual G4double GetP2Zmax(G4int id) const = 0;

    // Setters for attributes for plotting
    virtual G4bool SetP2Title(G4int id, const G4String& title) = 0;
    virtual G4bool SetP2XAxisTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetP2YAxisTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetP2ZAxisTitle(G4int id, const G4String& title) = 0;

    // Access attributes for plotting
    virtual G4String GetP2Title(G4int id) const = 0;
    virtual G4String GetP2XAxisTitle(G4int id) const = 0;
    virtual G4String GetP2YAxisTitle(G4int id) const = 0;
    virtual G4String GetP2ZAxisTitle(G4int id) const = 0;

    // Methods to manipulate histograms
    virtual G4bool WriteOnAscii(std::ofstream& output) = 0;

    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() = 0;
};

#endif

