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
// $Id: G4H1DummyManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for H1 for managers which do not support this object 
// type. It will just issue warnings. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4H1DummyManager_h
#define G4H1DummyManager_h 1

#include "G4VH1Manager.hh"
#include "globals.hh"

class G4H1DummyManager : public G4VH1Manager
{
  public:
    G4H1DummyManager(const G4AnalysisManagerState& state);
    virtual ~G4H1DummyManager();

  protected:
    // Virtual functions from base class
  
    // Methods to create histograms
    //
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binScheme = "linear");
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none");
                           
    virtual G4bool SetH1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binSchemeName = "linear");
    virtual G4bool SetH1(G4int id,
                           const std::vector<G4double>& edges,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none");

    virtual G4bool ScaleH1(G4int id, G4double factor);
                           
    // Method to fill histograms
    //
    virtual G4bool FillH1(G4int id, G4double value, G4double weight = 1.0);

    // Access methods
    // 
    virtual G4int  GetNofH1s() const;
    virtual G4int  GetH1Id(const G4String& name, G4bool warn = true) const;

    // Access to H1 parameters
    virtual G4int    GetH1Nbins(G4int id) const;
    virtual G4double GetH1Xmin(G4int id) const;
    virtual G4double GetH1Xmax(G4int id) const;
    virtual G4double GetH1Width(G4int id) const;

    // Attributes for plotting
    //

    // Setters
    virtual G4bool SetH1Title(G4int id, const G4String& title);
    virtual G4bool SetH1XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH1YAxisTitle(G4int id, const G4String& title);

    // Accessors
    virtual G4String GetH1Title(G4int id) const;
    virtual G4String GetH1XAxisTitle(G4int id) const;
    virtual G4String GetH1YAxisTitle(G4int id) const;

    // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output);
};

#endif

