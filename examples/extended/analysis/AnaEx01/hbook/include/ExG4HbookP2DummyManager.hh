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
// $Id: ExG4HbookP2DummyManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

/// \file hbook/include/ExG4HbookP2DummyManager.hh
/// \brief Definition of the ExG4HbookP2DummyManager class

// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef ExG4HbookP2DummyManager_h
#define ExG4HbookP2DummyManager_h 1

#include "G4VP2Manager.hh"

namespace tools {
namespace histo { 
class h3d; 
}
}

/// Manager class for P2 with dummy implementation. 
///
/// It will just issue warnings. 

class ExG4HbookP2DummyManager : public G4VP2Manager
{
  public:
    ExG4HbookP2DummyManager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookP2DummyManager();
   
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
                           const G4String& ybinScheme = "linear");
                           
    virtual G4int CreateP2(const G4String& name, const G4String& title,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none");
                          
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
                           const G4String& ybinScheme = "linear");
                           
    virtual G4bool SetP2(G4int id,
                           const std::vector<G4double>& xedges,
                           const std::vector<G4double>& yedges,
                           G4double zmin = 0, G4double zmax = 0,
                           const G4String& xunitName = "none", 
                           const G4String& yunitName = "none",
                           const G4String& zunitName = "none",
                           const G4String& xfcnName = "none", 
                           const G4String& yfcnName = "none",
                           const G4String& zfcnName = "none");

    virtual G4bool ScaleP2(G4int id, G4double factor);
    
    // Method to fill histograms
    //
    virtual G4bool FillP2(G4int id, 
                          G4double xvalue, G4double yvalue, G4double zvalue,
                          G4double weight = 1.0);
                          

    // Methods to manipulate histograms
    //

    // Access methods
    virtual G4int  GetP2Id(const G4String& name, G4bool warn = true) const;

    // Access to P2 parameters
    virtual G4int    GetP2Nxbins(G4int id) const;
    virtual G4double GetP2Xmin(G4int id) const;
    virtual G4double GetP2Xmax(G4int id) const;
    virtual G4double GetP2XWidth(G4int id) const;
    virtual G4int    GetP2Nybins(G4int id) const;
    virtual G4double GetP2Ymin(G4int id) const;
    virtual G4double GetP2Ymax(G4int id) const;
    virtual G4double GetP2YWidth(G4int id) const;
    virtual G4double GetP2Zmin(G4int id) const;
    virtual G4double GetP2Zmax(G4int id) const;
        
    // Setters for attributes for plotting
    virtual G4bool SetP2Title(G4int id, const G4String& title);
    virtual G4bool SetP2XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetP2YAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetP2ZAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    virtual G4String GetP2Title(G4int id) const;
    virtual G4String GetP2XAxisTitle(G4int id) const;
    virtual G4String GetP2YAxisTitle(G4int id) const;
    virtual G4String GetP2ZAxisTitle(G4int id) const;
 
     // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output);
    
  private:
    // methods
    void ExceptionForProfiles(const G4String& functionName); 
    void ExceptionForProfilesConst(const G4String& functionName) const; 

    // data members
    G4bool fWarn;    
};

#endif

