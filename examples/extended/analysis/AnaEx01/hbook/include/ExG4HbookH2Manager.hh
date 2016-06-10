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
//
/// \file hbook/include/ExG4HbookH2Manager.hh
/// \brief Definition of the ExG4HbookH2Manager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookH2Manager_h
#define ExG4HbookH2Manager_h 1

#include "G4VH2Manager.hh"
#include "ExG4HbookBaseHnManager.hh"
#include "globals.hh"

#include <tools/hbook/h2>

#include <vector>
#include <map>

class ExG4HbookFileManager;

struct h2_booking {
  h2_booking(G4int nxbins, G4double xmin, G4double xmax,
             G4int nybins, G4double ymin, G4double ymax)
    : fTitle(""),
      fNxbins(nxbins), 
      fXmin(xmin), 
      fXmax(xmax),
      fNybins(nybins), 
      fYmin(ymin), 
      fYmax(ymax) {}
  G4String fTitle;    
  G4int fNxbins;
  G4double fXmin;
  G4double fXmax;
  G4int fNybins;
  G4double fYmin;
  G4double fYmax;
};    

/// Manager class for HBook H2 histograms
///
/// The class implements the G4VH2Manager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookH2Manager : public G4VH2Manager
{
  friend class ExG4HbookAnalysisManager;

  protected:
     ExG4HbookH2Manager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookH2Manager();

    // Functions specific to the output type
    //

    // Set the offset of HBOOK ID for H2
    // ( default value = firstHistoID + 100 if firstHistoID > 0; 
    //   otherwise = 101 )
    G4bool SetH2HbookIdOffset(G4int offset);
    G4int  GetH2HbookIdOffset() const;
 
    // Set methods
    void SetFileManager(ExG4HbookFileManager* fileManager);

    // Access methods
    tools::hbook::h2*  GetH2(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
                              
    // Iterators
    std::vector<tools::hbook::h2*>::iterator BeginH2();
    std::vector<tools::hbook::h2*>::iterator EndH2();
    std::vector<tools::hbook::h2*>::const_iterator BeginConstH2() const;
    std::vector<tools::hbook::h2*>::const_iterator EndConstH2() const;
    
  protected:
    // Virtual functions from base class
    //

    // Methods to create histogrammes, ntuples
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
                           
    // Methods to fill histogrammes, ntuples
    virtual G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                          G4double weight = 1.0);
    
    // Access methods via names
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

  private:
    // static data members
    //
    static const G4int fgkDefaultH2HbookIdOffset;

    // methods
    //
    void SetH2HbookIdOffset();
    void CreateH2sFromBooking();

    void Reset();
    virtual h2_booking* GetH2Booking(G4int id, G4bool warn = true) const;

    virtual tools::hbook::h2*  GetH2InFunction(G4int id, G4String function,
                                      G4bool warn = true,
                                      G4bool onlyIfActive = true) const;

    // data members
    //
    ExG4HbookBaseHnManager fBaseToolsManager;
    ExG4HbookFileManager*  fFileManager;
    G4int fH2HbookIdOffset;
    std::vector<tools::hbook::h2*>  fH2Vector;            
    std::vector<h2_booking*>  fH2BookingVector;            
    std::map<G4String, G4int>  fH2NameIdMap;            
};

// inline functions

inline void ExG4HbookH2Manager::SetFileManager(ExG4HbookFileManager* fileManager)
{ fFileManager = fileManager; }

inline G4int ExG4HbookH2Manager::GetH2HbookIdOffset() const {
  return fH2HbookIdOffset;
}  

inline  std::vector<tools::hbook::h2*>::iterator ExG4HbookH2Manager::BeginH2()
{ return fH2Vector.begin(); }

inline  std::vector<tools::hbook::h2*>::iterator ExG4HbookH2Manager::EndH2()
{ return fH2Vector.end(); }

inline  std::vector<tools::hbook::h2*>::const_iterator 
ExG4HbookH2Manager::BeginConstH2() const
{ return fH2Vector.begin(); }

inline  std::vector<tools::hbook::h2*>::const_iterator 
ExG4HbookH2Manager::EndConstH2() const
{ return fH2Vector.end(); }


#endif 

#endif
