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
/// \file electromagnetic/TestEm4/include/ExG4HbookAnalysisManager.hh
/// \brief Definition of the ExG4HbookAnalysisManager class
//
// $Id$
//
/// \file ExG4HbookAnalysisManager.hh
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
#include <tools/hbook/wntuple>

#include <vector>
#include <map>

#define setpawc setpawc_
#define setntuc setntuc_
//#define ntuc ntuc_

class ExG4HbookAnalysisManager;

namespace G4Hbook {

  typedef tools::hbook::h1  G4AnaH1;
  typedef ExG4HbookAnalysisManager G4AnalysisManager; 
} 

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
    static ExG4HbookAnalysisManager* Instance();
  
    // Methods to manipulate files
    virtual G4bool OpenFile(const G4String& fileName);
    virtual G4bool Write();
    virtual G4bool CloseFile(); 
   
    // Methods to create histogrammes, ntuples
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax);
    virtual G4int CreateH2(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax);

    virtual void  CreateNtuple(const G4String& name, const G4String& title);
    virtual G4int CreateNtupleIColumn(const G4String& name);
    virtual G4int CreateNtupleFColumn(const G4String& name);
    virtual G4int CreateNtupleDColumn(const G4String& name);     
    virtual void  FinishNtuple();   
  
    // Methods to fill histogrammes, ntuples
    virtual G4bool FillH1(G4int id, G4double value, G4double weight = 1.0);
    virtual G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                          G4double weight = 1.0);
    virtual G4bool FillNtupleIColumn(G4int id, G4int value);
    virtual G4bool FillNtupleFColumn(G4int id, G4float value);
    virtual G4bool FillNtupleDColumn(G4int id, G4double value);
    virtual G4bool AddNtupleRow();
    
    // Access methods
    virtual tools::hbook::h1*  GetH1(G4int id, G4bool warn = true) const;
    virtual tools::hbook::h2*  GetH2(G4int id, G4bool warn = true) const;
    virtual tools::hbook::wntuple* GetNtuple() const;
    //tools::hbook::h1*  GetH1(const G4String& name, G4bool warn = true) const;
    
    // HBOOK does not allow IDs the same IDs for H1 and H2,
    // and also IDs starting from 0; thats why there is defined an offset
    // with respect to the G4AnalysisManager generic Ids.
    // The default values of these offsets can be changed by the user.
    //
    // Set the offset of HBOOK ID for H1
    // ( default value = firstHistoID if firstHistoID > 0; otherwise = 1)
    G4bool SetH1HbookIdOffset(G4int offset);
    //
    // Set the offset of HBOOK ID for H2
    // ( default value = firstHistoID + 100 if firstHistoID > 0; otherwise = 101 )
    G4bool SetH2HbookIdOffset(G4int offset);
    //
    // Set the HBOOK ID for the ntuple 
    // (default value = 1 )
    G4bool SetNtupleHbookId(G4int ntupleId);
    
    G4int  GetH1HbookIdOffset() const;
    G4int  GetH2HbookIdOffset() const;
    G4int  GetNtupleHbookId() const;
        
  private:
    // static data members
    //
    static ExG4HbookAnalysisManager* fgInstance;
    static const G4int fgkDefaultH2HbookIdOffset;
    static const G4int fgkDefaultNtupleHbookId;
    static const G4String fgkDefaultNtupleDirectoryName;
    
    // methods
    //
    tools::hbook::wntuple::column<int>*    GetNtupleIColumn(G4int id) const;
    tools::hbook::wntuple::column<float>*  GetNtupleFColumn(G4int id) const;
    tools::hbook::wntuple::column<double>* GetNtupleDColumn(G4int id) const;
 
    // data members
    //
    G4int fH1HbookIdOffset;
    G4int fH2HbookIdOffset;
    G4int fNtupleHbookId;

    tools::hbook::wfile*  fFile;
    
    std::vector<tools::hbook::h1*>         fH1Vector;            
    std::map<G4String, tools::hbook::h1*>  fH1MapByName;            

    std::vector<tools::hbook::h2*>         fH2Vector;            
    std::map<G4String, tools::hbook::h2*>  fH2MapByName;            

    G4String  fNtupleName;
    G4String  fNtupleTitle;
    tools::hbook::wntuple*  fNtuple; 
    std::map<G4int, tools::hbook::wntuple::column<int>* >    fNtupleIColumnMap;           
    std::map<G4int, tools::hbook::wntuple::column<float>* >  fNtupleFColumnMap;           
    std::map<G4int, tools::hbook::wntuple::column<double>* > fNtupleDColumnMap;           
};

// inline functions

inline G4int ExG4HbookAnalysisManager::GetH1HbookIdOffset() const {
  return fH1HbookIdOffset;
}  

inline G4int ExG4HbookAnalysisManager::GetH2HbookIdOffset() const {
  return fH2HbookIdOffset;
}  

inline G4int ExG4HbookAnalysisManager::GetNtupleHbookId() const {
  return fNtupleHbookId;
}  

#endif 

#endif
