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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifndef G4RootAnalysisManager_h
#define G4RootAnalysisManager_h 1

#include "G4VAnalysisManager.hh"
#include "globals.hh"

#include "tools/wroot/file"
#include "tools/wroot/to"
#include "tools/wroot/ntuple"
#include "tools/histo/h1d"
#include "tools/histo/h2d"

#include <vector>
#include <map>

class G4RootAnalysisManager : public G4VAnalysisManager
{
  public:
    G4RootAnalysisManager();
    virtual ~G4RootAnalysisManager();
    
    // static methods
    static G4RootAnalysisManager* Instance();

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
    virtual tools::histo::h1d*  GetH1(G4int id, G4bool warn = true) const;
    virtual tools::histo::h2d*  GetH2(G4int id, G4bool warn = true) const;
    //tools::histo::h1d*  GetH1D(const G4String& name, G4bool warn = true) const;
        // later
        
  private:
    // static data members
    //
    static G4RootAnalysisManager* fgInstance;

    // methods
    //
    tools::wroot::ntuple::column<int>*    GetNtupleIColumn(G4int id) const;
    tools::wroot::ntuple::column<float>*  GetNtupleFColumn(G4int id) const;
    tools::wroot::ntuple::column<double>* GetNtupleDColumn(G4int id) const;
 
    // data members
    //
    tools::wroot::file*       fFile;
    tools::wroot::directory*  fHistoDirectory;
    tools::wroot::directory*  fNtupleDirectory;

    std::vector<tools::histo::h1d*>         fH1Vector;            
    std::map<G4String, tools::histo::h1d*>  fH1MapByName;            
    
    std::vector<tools::histo::h2d*>         fH2Vector;            
    std::map<G4String, tools::histo::h2d*>  fH2MapByName;            
    
    G4String  fNtupleName;
    G4String  fNtupleTitle;
    tools::wroot::ntuple*  fNtuple; 
    std::map<G4int, tools::wroot::ntuple::column<int>* >    fNtupleIColumnMap;           
    std::map<G4int, tools::wroot::ntuple::column<float>* >  fNtupleFColumnMap;           
    std::map<G4int, tools::wroot::ntuple::column<double>* > fNtupleDColumnMap;           
};

#endif

