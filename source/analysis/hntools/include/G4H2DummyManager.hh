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
// $Id: G4H2DummyManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for H2 for managers which do not support this object 
// type. It will just issue warnings. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4H2DummyManager_h
#define G4H2DummyManager_h 1

#include "G4VH2Manager.hh"
#include "globals.hh"

class G4H2DummyManager : public G4VH2Manager
{
  public:
    G4H2DummyManager(const G4AnalysisManagerState& state);
    virtual ~G4H2DummyManager();

  protected:
    // Virtual functions from base class
    //

    // Methods to create histograms
    //
    virtual G4int CreateH2(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear");
                           
    virtual G4int CreateH2(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none");
                          
    virtual G4bool SetH2(G4int id,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& xbinScheme = "linear",
                           const G4String& ybinScheme = "linear");
                           
    virtual G4bool SetH2(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none");

    virtual G4bool ScaleH2(G4int id, G4double factor);

    // Method to fill histograms
    //
    virtual G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                          G4double weight = 1.0);
                          
    // Methods to manipulate histograms
    //

    // Access methods
    virtual G4int  GetNofH2s() const;
    virtual G4int  GetH2Id(const G4String& name, G4bool warn = true) const;

    // Access to H2 parameters
    virtual G4int    GetH2Nxbins(G4int id) const;
    virtual G4double GetH2Xmin(G4int id) const;
    virtual G4double GetH2Xmax(G4int id) const;
    virtual G4double GetH2XWidth(G4int id) const;
    virtual G4int    GetH2Nybins(G4int id) const;
    virtual G4double GetH2Ymin(G4int id) const;
    virtual G4double GetH2Ymax(G4int id) const;
    virtual G4double GetH2YWidth(G4int id) const;
        
    // Setters for attributes for plotting
    virtual G4bool SetH2Title(G4int id, const G4String& title);
    virtual G4bool SetH2XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH2YAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH2ZAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    virtual G4String GetH2Title(G4int id) const;
    virtual G4String GetH2XAxisTitle(G4int id) const;
    virtual G4String GetH2YAxisTitle(G4int id) const;
    virtual G4String GetH2ZAxisTitle(G4int id) const;
 
     // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output);
};

#endif

