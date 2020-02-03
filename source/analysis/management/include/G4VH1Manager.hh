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

// Base class for H1 manager.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VH1Manager_h
#define G4VH1Manager_h 1

#include "globals.hh"

#include <vector>
#include <memory>

class G4HnManager;

class G4VH1Manager
{
  // Disable using the object managers outside 
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;

  public:
    G4VH1Manager() {}
    virtual ~G4VH1Manager() {}

    // deleted copy constructor & assignment operator
    G4VH1Manager(const G4VH1Manager& rhs) = delete;
    G4VH1Manager& operator=(const G4VH1Manager& rhs) = delete;

  protected:
    // Methods for handling histograms
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binSchemeName = "linear") = 0;
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none") = 0;

    virtual G4bool SetH1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binSchemeName = "linear") = 0;
    virtual G4bool SetH1(G4int id,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none") = 0;

    virtual G4bool ScaleH1(G4int id, G4double factor) = 0;
                           
    // Methods to fill histograms
    virtual G4bool FillH1(G4int id, G4double value, G4double weight = 1.0) = 0;

    // Access methods
    virtual G4int  GetH1Id(const G4String& name, G4bool warn = true) const = 0;

    // Access to H1 parameters
    virtual G4int    GetH1Nbins(G4int id) const = 0;
    virtual G4double GetH1Xmin(G4int id) const = 0;
    virtual G4double GetH1Xmax(G4int id) const = 0;
    virtual G4double GetH1Width(G4int id) const = 0;

    // Setters for attributes for plotting
    virtual G4bool SetH1Title(G4int id, const G4String& title) = 0;
    virtual G4bool SetH1XAxisTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetH1YAxisTitle(G4int id, const G4String& title) = 0;

    // Access attributes for plotting
    virtual G4String GetH1Title(G4int id) const = 0;
    virtual G4String GetH1XAxisTitle(G4int id) const = 0;
    virtual G4String GetH1YAxisTitle(G4int id) const = 0;

    // Methods to manipulate histograms
    virtual G4bool WriteOnAscii(std::ofstream& output) = 0;

    // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() = 0;
};

#endif

