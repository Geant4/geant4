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
//
// Author: Ivana Hrivnacova, 30/10/2018  (ivana@ipno.in2p3.fr)

#ifndef G4ScoreNtupleWriter_h
#define G4ScoreNtupleWriter_h 1

#include "G4VScoreNtupleWriter.hh"
#include "G4HCofThisEvent.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Threading.hh"
#include "globals.hh"

#include <memory>

class G4HCofThisEvent;
class G4VAnalysisManager;
class G4ScoreNtupleWriterMessenger;

// class description:
//
// This class implements storing hits collections of G4THitsMap<G4double>
// type vith Geant4 analysis tools.
//
// An ntuple with three columns is created for each primitive scorer:
// int column - eventNumber 
// int column - copyNumber
// double column - scored value

class G4ScoreNtupleWriter : public G4VScoreNtupleWriter
{
  public: 
    G4ScoreNtupleWriter(const G4String& outputType = "none");
    virtual ~G4ScoreNtupleWriter();
    
    // methods
    virtual G4bool Book(G4HCofThisEvent* hce);
    virtual void OpenFile();
    virtual void Fill(G4HCofThisEvent* hce, G4int eventNumber);
    virtual void Write();

    // set methods
    void SetOutputType(const G4String& outputTypeName);
    void SetFileName(const G4String& fileName);
    void SetVerboseLevel(G4int value);

    // get methods
    G4AnalysisOutput GetOutputType() const;
    G4String GetFileName() const   { return fFileName; }
    G4int GetVerboseLevel() const  { return fVerboseLevel; }

  protected:
    // methods
    virtual G4VScoreNtupleWriter* CreateInstance() const;

  private:
    // methods
    void CreateAnalysisManager();

    // data members
    G4ScoreNtupleWriterMessenger* fMessenger;
    std::vector<G4int>  fHCIds;
    G4VAnalysisManager* fAnalysisManager;
    G4String fOutputTypeName;
    G4String fFileName;
    G4int    fVerboseLevel;
    G4bool   fHasAnalysisManager;
    G4bool   fHasAnalysisFile;
    G4bool   fIsBooked;
    G4bool   fIsInitialized;
    G4int    fFirstNtupleId;
};

#endif
