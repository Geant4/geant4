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
/// \file hbook/include/ExG4HbookAnalysisManager.hh
/// \brief Definition of the ExG4HbookAnalysisManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookAnalysisManager_h
#define ExG4HbookAnalysisManager_h 1

#include "G4VAnalysisManager.hh"
#include "globals.hh"

#include <tools/hbook/wfile>
#include <tools/hbook/h1>
#include <tools/hbook/h2>
#include <tools/hbook/p1>
#include <tools/hbook/wntuple>

#include <vector>
#include <map>

#define setpawc setpawc_
#define setntuc setntuc_
//#define ntuc ntuc_

class ExG4HbookAnalysisManager;

namespace G4Hbook {
  typedef tools::hbook::h1  G4AnaH1;
  typedef tools::hbook::h1  G4H1;
  typedef tools::hbook::h2  G4AnaH2;
  typedef tools::hbook::h2  G4H2;
  typedef tools::hbook::p1  G4P1;
  typedef tools::hbook::wntuple  G4Ntuple; 
  typedef ExG4HbookAnalysisManager G4AnalysisManager; 
} 

class ExG4HbookFileManager;
class ExG4HbookH1Manager;
class ExG4HbookH2Manager;
class ExG4HbookH3DummyManager;
class ExG4HbookP1Manager;
class ExG4HbookP2DummyManager;
class ExG4HbookNtupleManager;

/// HBook Analysis manager
///
/// The class implements the G4VAnalysisManager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookAnalysisManager : public G4VAnalysisManager
{
  public:
    ExG4HbookAnalysisManager();
    virtual ~ExG4HbookAnalysisManager();
    
    // static methods
    static ExG4HbookAnalysisManager* Create(G4bool isMaster = true);
    static ExG4HbookAnalysisManager* Instance();
  
    // HBOOK does not allow IDs the same IDs for H1 and H2,
    // and also IDs starting from 0; thats why there is defined an offset
    // with respect to the G4AnalysisManager generic Ids.
    // The default values of these offsets can be changed by the user.
    //
    // Set the offset of HBOOK ID for H1
    // ( default value = firstHistoID if firstHistoID > 0; otherwise = 1)
    G4bool SetH1HbookIdOffset(G4int offset);
    G4int  GetH1HbookIdOffset() const;

    // Set the offset of HBOOK ID for H2
    // ( default value = firstHistoID + 100 if firstHistoID > 0; 
    //   otherwise = 101 )
    G4bool SetH2HbookIdOffset(G4int offset);
    G4int  GetH2HbookIdOffset() const;

    // Set the offset of NTUPLE ID for ntuples
    // ( default value = firstNtupleID if firstNtupleID > 0; otherwise = 1)
    G4bool SetNtupleHbookIdOffset(G4int offset);
    G4int  GetNtupleHbookIdOffset() const;

    // Access methods
    //
    tools::hbook::h1*  GetH1(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::hbook::h2*  GetH2(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::hbook::wntuple* GetNtuple() const;
    tools::hbook::wntuple* GetNtuple(G4int ntupleId) const;
 
     // Iterators
    std::vector<tools::hbook::h1*>::iterator BeginH1();
    std::vector<tools::hbook::h1*>::iterator EndH1();
    std::vector<tools::hbook::h1*>::const_iterator BeginConstH1() const;
    std::vector<tools::hbook::h1*>::const_iterator EndConstH1() const;
    
    std::vector<tools::hbook::h2*>::iterator BeginH2();
    std::vector<tools::hbook::h2*>::iterator EndH2();
    std::vector<tools::hbook::h2*>::const_iterator BeginConstH2() const;
    std::vector<tools::hbook::h2*>::const_iterator EndConstH2() const;
    
    std::vector<tools::hbook::wntuple*>::iterator BeginNtuple();
    std::vector<tools::hbook::wntuple*>::iterator EndNtuple();
    std::vector<tools::hbook::wntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::hbook::wntuple*>::const_iterator EndConstNtuple() const;
    
  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName);
    virtual G4bool WriteImpl();
    virtual G4bool CloseFileImpl(); 
   
  private:
    // static data members
    //
    static ExG4HbookAnalysisManager* fgInstance;
 
    // methods
    //
    void Reset();
   
    ExG4HbookH1Manager*      fH1Manager;
    ExG4HbookH2Manager*      fH2Manager;
    ExG4HbookH3DummyManager* fH3Manager;
    ExG4HbookP1Manager*      fP1Manager;
    ExG4HbookP2DummyManager* fP2Manager;
    ExG4HbookNtupleManager*  fNtupleManager;
    ExG4HbookFileManager*    fFileManager;
};

#include "ExG4HbookAnalysisManager.icc"

#endif 

#endif
