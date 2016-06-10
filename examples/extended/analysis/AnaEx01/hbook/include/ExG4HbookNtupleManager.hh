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
/// \file hbook/include/ExG4HbookNtupleManager.hh
/// \brief Definition of the ExG4HbookNtupleManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookNtupleManager_h
#define ExG4HbookNtupleManager_h 1

#include "G4VNtupleManager.hh"
#include "ExG4HbookNtupleDescription.hh"
#include "globals.hh"

#include <tools/hbook/wntuple>

#include <vector>
#include <map>

class ExG4HbookFileManager;

/// Manager class for HBook ntuples
///
/// The class implements the G4VAnalysisManager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookNtupleManager : public G4VNtupleManager
{
  friend class ExG4HbookAnalysisManager;

  protected:
    ExG4HbookNtupleManager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookNtupleManager();

    // Functions specific to the output type
    //

    // Set the offset of NTUPLE ID for ntuples
    // ( default value = firstNtupleID if firstNtupleID > 0; otherwise = 1)
    G4bool SetNtupleHbookIdOffset(G4int offset);
    G4int  GetNtupleHbookIdOffset() const;

    // Set methods
    void SetFileManager(ExG4HbookFileManager* fileManager);

    // Access methods
    tools::hbook::wntuple* GetNtuple() const;
    tools::hbook::wntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::hbook::wntuple*>::iterator BeginNtuple();
    std::vector<tools::hbook::wntuple*>::iterator EndNtuple();
    std::vector<tools::hbook::wntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::hbook::wntuple*>::const_iterator EndConstNtuple() const;

    // Virtual functions from base class
    //
    // Methods to create ntuples
    virtual G4int CreateNtuple(const G4String& name, const G4String& title);
    // Create columns in the last created ntuple
    virtual G4int CreateNtupleIColumn(
                    const G4String& name, std::vector<int>* vector);
    virtual G4int CreateNtupleFColumn(
                    const G4String& name, std::vector<float>* vector);
    virtual G4int CreateNtupleDColumn(
                    const G4String& name, std::vector<double>* vector);
    virtual G4int CreateNtupleSColumn(const G4String& name);
    virtual void  FinishNtuple();   
    // Create columns in the ntuple with given id
    virtual G4int CreateNtupleIColumn(G4int ntupleId, 
                    const G4String& name, std::vector<int>* vector);
    virtual G4int CreateNtupleFColumn(G4int ntupleId, 
                    const G4String& name, std::vector<float>* vector);
    virtual G4int CreateNtupleDColumn(G4int ntupleId, 
                    const G4String& name, std::vector<double>* vector);
    virtual G4int CreateNtupleSColumn(G4int ntupleId, const G4String& name);
    virtual void  FinishNtuple(G4int ntupleId);   
 
    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId                     
    virtual G4bool FillNtupleIColumn(G4int columnId, G4int value);
    virtual G4bool FillNtupleFColumn(G4int columnId, G4float value);
    virtual G4bool FillNtupleDColumn(G4int columnId, G4double value);
    virtual G4bool FillNtupleSColumn(G4int columnId, const G4String& value);
    virtual G4bool AddNtupleRow();
    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)
    virtual G4bool FillNtupleIColumn(
                       G4int ntupleId, G4int columnId, G4int value);
    virtual G4bool FillNtupleFColumn(
                       G4int ntupleId, G4int columnId, G4float value);
    virtual G4bool FillNtupleDColumn(
                       G4int ntupleId, G4int columnId, G4double value);
    virtual G4bool FillNtupleSColumn(
                       G4int ntupleId, G4int columnId, const G4String& value);
    virtual G4bool AddNtupleRow(G4int ntupleId);
    
    // Access methods
    virtual G4int GetNofNtuples() const;

  private:
    // methods
    //
    void SetNtupleHbookIdOffset();
    void CreateNtuplesFromBooking();

    tools::hbook::wntuple::column<int>*    
      GetNtupleIColumn(G4int ntupleId, G4int columnId) const;
    tools::hbook::wntuple::column<float>*  
      GetNtupleFColumn(G4int ntupleId, G4int columnId) const;
    tools::hbook::wntuple::column<double>* 
      GetNtupleDColumn(G4int ntupleId, G4int columnId) const;

    void Reset();
    ExG4HbookNtupleDescription*  GetNtupleInFunction(G4int id, 
                                           G4String function,
                                           G4bool warn = true,
                                           G4bool onlyIfActive = true) const;
                                           
    // data members
    //
    ExG4HbookFileManager* fFileManager;
    G4int fNtupleHbookIdOffset;
    std::vector<ExG4HbookNtupleDescription*> fNtupleDescriptionVector;
    std::vector<tools::hbook::wntuple*> fNtupleVector;
};

// inline functions

inline void ExG4HbookNtupleManager::SetFileManager(ExG4HbookFileManager* fileManager)
{ fFileManager = fileManager; }

inline G4int ExG4HbookNtupleManager::GetNtupleHbookIdOffset() const {
  return fNtupleHbookIdOffset;
}  

inline std::vector<tools::hbook::wntuple*>::iterator 
ExG4HbookNtupleManager::BeginNtuple()
{ return fNtupleVector.begin(); }

inline std::vector<tools::hbook::wntuple*>::iterator 
ExG4HbookNtupleManager::EndNtuple()
{ return fNtupleVector.end(); }

inline std::vector<tools::hbook::wntuple*>::const_iterator 
ExG4HbookNtupleManager::BeginConstNtuple() const
{ return fNtupleVector.begin(); }

inline  std::vector<tools::hbook::wntuple*>::const_iterator 
ExG4HbookNtupleManager::EndConstNtuple() const
{ return fNtupleVector.end(); }

inline G4int ExG4HbookNtupleManager::GetNofNtuples() const
{ return fNtupleVector.size(); }


#endif 

#endif
