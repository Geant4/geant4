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
/// \file hbook/include/ExG4HbookH1Manager.hh
/// \brief Definition of the ExG4HbookH1Manager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookH1Manager_h
#define ExG4HbookH1Manager_h 1

#include "G4VH1Manager.hh"
#include "G4HnInformation.hh"
#include "ExG4HbookBaseHnManager.hh"
#include "globals.hh"

#include <tools/hbook/h1>

#include <vector>
#include <map>

class ExG4HbookFileManager;

struct h1_booking {
  h1_booking(G4int nbins, G4double xmin, G4double xmax)
    : fTitle(""),
      fNbins(nbins), 
      fXmin(xmin), 
      fXmax(xmax),
      fEdges() {}
  h1_booking(const std::vector<G4double>& edges)
    : fTitle(""),
      fNbins(0), 
      fXmin(0), 
      fXmax(0),
      fEdges() {
    for (G4int i=0; i<=G4int(edges.size()); ++i) fEdges.push_back(edges[i]);
  }
  G4String fTitle;    
  G4int fNbins;
  G4double fXmin;
  G4double fXmax;
  std::vector<G4double> fEdges;
};  
  
/// Manager class for HBook H1 histograms
///
/// The class implements the G4VH1Manager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookH1Manager : public G4VH1Manager
{
  friend class ExG4HbookAnalysisManager;

  protected:
    ExG4HbookH1Manager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookH1Manager();

    // Functions specific to the output type
    //

    // HBOOK does not allow IDs the same IDs for H1 and H2,
    // and also IDs starting from 0; thats why there is defined an offset
    // with respect to the G4AnalysisManager generic Ids.
    // The default values of these offsets can be changed by the user.
    //
    // Set the offset of HBOOK ID for H1
    // ( default value = firstHistoID if firstHistoID > 0; otherwise = 1)
    G4bool SetH1HbookIdOffset(G4int offset);
    G4int  GetH1HbookIdOffset() const;
 
    // Set methods
    void SetFileManager(ExG4HbookFileManager* fileManager);

    // Access methods
    //
    tools::hbook::h1*  GetH1(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
                              
    // Iterators
    std::vector<tools::hbook::h1*>::iterator BeginH1();
    std::vector<tools::hbook::h1*>::iterator EndH1();
    std::vector<tools::hbook::h1*>::const_iterator BeginConstH1() const;
    std::vector<tools::hbook::h1*>::const_iterator EndConstH1() const;
    
    // Virtual functions from base class
    //

    // Methods to create histogrammes, ntuples
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           const G4String& unitName = "none",
                           const G4String& fcnName = "none",
                           const G4String& binSchemeName = "linear");
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
                            
                           
    // Methods to fill histogrammes, ntuples
    virtual G4bool FillH1(G4int id, G4double value, G4double weight = 1.0);

    // Access methods
    //
    virtual G4int  GetH1Id(const G4String& name, G4bool warn = true) const;

   // Access to H1 parameters
    virtual G4int    GetH1Nbins(G4int id) const;
    virtual G4double GetH1Xmin(G4int id) const;
    virtual G4double GetH1Xmax(G4int id) const;
    virtual G4double GetH1Width(G4int id) const;
        
    // Setters for attributes for plotting
    virtual G4bool SetH1Title(G4int id, const G4String& title);
    virtual G4bool SetH1XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetH1YAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    virtual G4String GetH1Title(G4int id) const;
    virtual G4String GetH1XAxisTitle(G4int id) const;
    virtual G4String GetH1YAxisTitle(G4int id) const;

    // Write data on ASCII file
    virtual G4bool WriteOnAscii(std::ofstream& output);

  private:
    // methods
    //
    void SetH1HbookIdOffset();
    void AddH1Information(const G4String& name,  
                          const G4String& unitName, 
                          const G4String& fcnName,
                          G4BinScheme binScheme) const;
    
    G4int CreateH1FromBooking(h1_booking* h1Booking, 
                          G4bool chDir = true);
    G4int RegisterH1Booking(const G4String& name, 
                          h1_booking* h1Booking);

    void  BeginCreateH1(const G4String& name);
    G4int FinishCreateH1(const G4String& name, h1_booking* h1Booking,
                         const G4String& unitName, const G4String& fcnName,
                         G4BinScheme binScheme);
    
    G4bool BeginSetH1(G4int id,
                         h1_booking* h1Booking,
                         G4HnInformation* info);
    G4bool FinishSetH1(G4int id,
                         G4HnInformation* info,
                         const G4String& unitName, const G4String& fcnName,
                         G4BinScheme binScheme);

    void CreateH1sFromBooking();
    void Reset();
    virtual h1_booking* GetH1Booking(G4int id, G4bool warn = true) const;

    virtual tools::hbook::h1*  GetH1InFunction(G4int id, G4String function,
                                      G4bool warn = true,
                                      G4bool onlyIfActive = true) const;

    // data members
    //
    ExG4HbookBaseHnManager fBaseToolsManager;
    ExG4HbookFileManager*  fFileManager;
    G4int fH1HbookIdOffset;
    std::vector<tools::hbook::h1*>  fH1Vector;            
    std::vector<h1_booking*>  fH1BookingVector;            
    std::map<G4String, G4int>  fH1NameIdMap;            
};

// inline functions

inline void ExG4HbookH1Manager::SetFileManager(ExG4HbookFileManager* fileManager)
{ fFileManager = fileManager; }

inline G4int ExG4HbookH1Manager::GetH1HbookIdOffset() const {
  return fH1HbookIdOffset;
}  

inline  std::vector<tools::hbook::h1*>::iterator ExG4HbookH1Manager::BeginH1()
{ return fH1Vector.begin(); }

inline  std::vector<tools::hbook::h1*>::iterator ExG4HbookH1Manager::EndH1()
{ return fH1Vector.end(); }

inline  std::vector<tools::hbook::h1*>::const_iterator 
ExG4HbookH1Manager::BeginConstH1() const
{ return fH1Vector.begin(); }

inline  std::vector<tools::hbook::h1*>::const_iterator 
ExG4HbookH1Manager::EndConstH1() const
{ return fH1Vector.end(); }


#endif 

#endif
