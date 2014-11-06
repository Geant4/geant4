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
/// \file hbook/include/ExG4HbookP1Manager.hh
/// \brief Definition of the ExG4HbookP1Manager class

// Author: Ivana Hrivnacova, 03/11/2014  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookP1Manager_h
#define ExG4HbookP1Manager_h 1

#include "G4VP1Manager.hh"
#include "G4HnInformation.hh"
#include "ExG4HbookBaseHnManager.hh"
#include "globals.hh"

#include <tools/hbook/p1>

#include <vector>
#include <map>

class ExG4HbookFileManager;

struct p1_booking {
  p1_booking(G4int nbins, G4double xmin, G4double xmax, 
             G4double ymin, G4double ymax)
    : fTitle(""),
      fNbins(nbins), 
      fXmin(xmin), 
      fXmax(xmax),
      fYmin(ymin), 
      fYmax(ymax),
      fEdges() {}
  p1_booking(const std::vector<G4double>& edges,
             G4double ymin, G4double ymax)
    : fTitle(""),
      fNbins(0), 
      fXmin(0), 
      fXmax(0),
      fYmin(ymin), 
      fYmax(ymax),
      fEdges() {
    for (G4int i=0; i<=G4int(edges.size()); ++i) fEdges.push_back(edges[i]);
  }
  G4String fTitle;    
  G4int fNbins;
  G4double fXmin;
  G4double fXmax;
  G4double fYmin;
  G4double fYmax;
  std::vector<G4double> fEdges;
};  
  
/// Manager class for HBook P1 histograms
///
/// The class implements the G4VP1Manager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookP1Manager : public G4VP1Manager
{
  friend class ExG4HbookAnalysisManager;

  protected:
    ExG4HbookP1Manager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookP1Manager();

    // Functions specific to the output type
    //

    // HBOOK does not allow IDs the same IDs for P1 and H2,
    // and also IDs starting from 0; thats why there is defined an offset
    // with respect to the G4AnalysisManager generic Ids.
    // The default values of these offsets can be changed by the user.
    //
    // Set the offset of HBOOK ID for P1
    // ( default value = firstHistoID if firstHistoID > 0; otherwise = 1)
    G4bool SetP1HbookIdOffset(G4int offset);
    G4int  GetP1HbookIdOffset() const;
 
    // Set methods
    void SetFileManager(ExG4HbookFileManager* fileManager);

    // Access methods
    //
    tools::hbook::p1*  GetP1(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
                              
    // Iterators
    std::vector<tools::hbook::p1*>::iterator BeginP1();
    std::vector<tools::hbook::p1*>::iterator EndP1();
    std::vector<tools::hbook::p1*>::const_iterator BeginConstP1() const;
    std::vector<tools::hbook::p1*>::const_iterator EndConstP1() const;
    
    // Virtual functions from base class
    //

    // Methods to create profiles
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
                            
                           
    // Methods to fill profiles
    virtual G4bool FillP1(G4int id, G4double value, G4double yvalue, 
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
        
    // Setters for attributes for plotting
    virtual G4bool SetP1Title(G4int id, const G4String& title);
    virtual G4bool SetP1XAxisTitle(G4int id, const G4String& title);
    virtual G4bool SetP1YAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    virtual G4String GetP1Title(G4int id) const;
    virtual G4String GetP1XAxisTitle(G4int id) const;
    virtual G4String GetP1YAxisTitle(G4int id) const;

    // Write data on ASCII file
    // virtual G4bool WriteOnAscii(std::ofstream& output);

  private:
    // methods
    //
    void SetP1HbookIdOffset();
    void AddP1Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme) const;
    
    G4int CreateP1FromBooking(p1_booking* p1Booking, 
                          G4bool chDir = true);
    G4int RegisterP1Booking(const G4String& name, 
                          p1_booking* p1Booking);

    void  BeginCreateP1(const G4String& name);
    G4int FinishCreateP1(const G4String& name, p1_booking* p1Booking,
                         const G4String& xunitName, 
                         const G4String& yunitName, 
                         const G4String& xfcnName,
                         const G4String& yfcnName,
                         G4BinScheme xbinScheme) ;
    
    G4bool BeginSetP1(G4int id,
                         p1_booking* p1Booking,
                         G4HnInformation* info);
    G4bool FinishSetP1(G4int id,
                         G4HnInformation* info,
                         const G4String& xunitName, 
                         const G4String& yunitName, 
                         const G4String& xfcnName,
                         const G4String& yfcnName,
                         G4BinScheme xbinScheme) ;

    void CreateP1sFromBooking();
    void Reset();
    virtual p1_booking* GetP1Booking(G4int id, G4bool warn = true) const;

    virtual tools::hbook::p1*  GetP1InFunction(G4int id, G4String function,
                                      G4bool warn = true,
                                      G4bool onlyIfActive = true) const;

    // data members
    //
    ExG4HbookBaseHnManager fBaseToolsManager;
    ExG4HbookFileManager*  fFileManager;
    G4int fP1HbookIdOffset;
    std::vector<tools::hbook::p1*>  fP1Vector;            
    std::vector<p1_booking*>  fP1BookingVector;            
    std::map<G4String, G4int>  fP1NameIdMap;            
};

// inline functions

inline void ExG4HbookP1Manager::SetFileManager(ExG4HbookFileManager* fileManager)
{ fFileManager = fileManager; }

inline G4int ExG4HbookP1Manager::GetP1HbookIdOffset() const 
{ return fP1HbookIdOffset; }  

inline  std::vector<tools::hbook::p1*>::iterator ExG4HbookP1Manager::BeginP1()
{ return fP1Vector.begin(); }

inline  std::vector<tools::hbook::p1*>::iterator ExG4HbookP1Manager::EndP1()
{ return fP1Vector.end(); }

inline  std::vector<tools::hbook::p1*>::const_iterator 
ExG4HbookP1Manager::BeginConstP1() const
{ return fP1Vector.begin(); }

inline  std::vector<tools::hbook::p1*>::const_iterator 
ExG4HbookP1Manager::EndConstP1() const
{ return fP1Vector.end(); }


#endif 

#endif
