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

// Manager class for P1 for managers which do not support this object 
// type. It will just issue warnings. 
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4P1DummyManager_h
#define G4P1DummyManager_h 1

#include "G4VP1Manager.hh"
#include "G4HnManager.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include <vector>
#include <map>

namespace tools {
namespace histo { 
class p1d; 
}
}

class G4P1DummyManager : public G4VP1Manager
{
  public:
    G4P1DummyManager(const G4AnalysisManagerState& state);
    virtual ~G4P1DummyManager();

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
                           const G4String& xbinScheme = "linear");
    virtual G4int CreateP1(const G4String& name, const G4String& title,
                           const std::vector<G4double>& edges,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none");
                           
    virtual G4bool SetP1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& xbinScheme = "linear");
    virtual G4bool SetP1(G4int id,
                           const std::vector<G4double>& edges,
                           G4double ymin = 0, G4double ymax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none");

    virtual G4bool ScaleP1(G4int id, G4double factor);
    
    // Method to fill profiles
    //
    virtual G4bool FillP1(G4int id, G4double xvalue, G4double yvalue, 
                          G4double weight = 1.0);

    // Access methods
    //
    virtual G4int  GetP1Id(const G4String& name, G4bool warn = true) const;

    // Access to P1 parameters
    virtual G4int    GetP1Nbins(G4int id) const;
    virtual G4double GetP1Xmin(G4int id) const;
    virtual G4double GetP1Xmax(G4int id) const;
    virtual G4double GetP1XWidth(G4int id) const;
    virtual G4double GetP1Ymin(G4int id) const;
    virtual G4double GetP1Ymax(G4int id) const;

    // Attributes for plotting
    //

    // Setters
    virtual G4bool SetP1Title(G4int id, const G4String& title);
    virtual G4bool SetP1XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetP1YAxisTitle(G4int id, const G4String& title);

    // Accessors
    virtual G4String GetP1Title(G4int id) const;
    virtual G4String GetP1XAxisTitle(G4int id) const;
    virtual G4String GetP1YAxisTitle(G4int id) const;

    // Write data on ASCII file
    //virtual G4bool WriteOnAscii(std::ofstream& output);
};

#endif

