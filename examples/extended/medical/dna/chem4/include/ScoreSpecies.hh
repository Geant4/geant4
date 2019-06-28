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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
#ifndef CHEM4_ScoreSpecies_h
#define CHEM4_ScoreSpecies_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include <set>

class G4VAnalysisManager;
class G4MolecularConfiguration;

//#define _ScoreSpecies_FOR_ALL_EVENTS

/** \file ScoreSpecies.hh*/

// Description:
//   This is a primitive scorer class for scoring the radiolitic species
// produced after irradiation in a water volume
//
// Created: 2015-10-27  by M. Karamitros,
// modified: 2016-03-16 by P. Piersimoni

class ScoreSpecies : public G4VPrimitiveScorer
{
public:
  ScoreSpecies(G4String name, G4int depth=0);

  virtual ~ScoreSpecies();
  
  /** Add a time at which the number of species should be recorded.
      Default times are set up to 1 microsecond.*/
  inline void AddTimeToRecord(double time)
  {
    fTimeToRecord.insert(time);
  }
  
  /**  Remove all times to record, must be reset by user.*/
  inline void ClearTimeToRecord()
  {
    fTimeToRecord.clear();
  }

  /** Get number of recorded events*/
  inline int GetNumberOfRecordedEvents() const
  {
    return fNEvent;
  }
  
  /** Write results to an text file*/
  void ASCII();
  
  /** Write results to whatever chosen file format*/
  void WriteWithAnalysisManager(G4VAnalysisManager*);
  
  struct SpeciesInfo
  {
    SpeciesInfo()
    {
      fNumber = 0;
      fG = 0.;
      fG2 = 0.;
    }
    SpeciesInfo(const SpeciesInfo& right) // Species A(B);
    {
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
    }
    SpeciesInfo& operator=(const SpeciesInfo& right) // A = B
    {
      if(&right == this) return *this;
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
      return *this;
    }
    int fNumber;
    double fG;
    double fG2;
  };
  
#ifdef _ScoreSpecies_FOR_ALL_EVENTS
  struct SpeciesInfoSOA
  {
    SpeciesInfoSOA()
    {
      //fNumber = 0;
      //fG = 0.;
      //fG2 = 0.;
    }
    SpeciesInfoSOA(const SpeciesInfoSOA& right) // Species A(B);
    {
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
      fEventID = right.fEventID;
    }
    SpeciesInfoSOA& operator=(const SpeciesInfoSOA& right) // A = B
    {
      if(&right == this) return *this;
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
      fEventID = right.fEventID;
      return *this;
    }
    std::vector<int> fNumber;
    std::vector<double> fG;
    std::vector<double> fG2;
    std::vector<int> fEventID;
  };
#endif
  
private:
  typedef const G4MolecularConfiguration Species;
  typedef std::map<Species*, SpeciesInfo>  InnerSpeciesMap;
  typedef std::map<double, InnerSpeciesMap> SpeciesMap;
  SpeciesMap fSpeciesInfoPerTime;

#ifdef _ScoreSpecies_FOR_ALL_EVENTS
  typedef std::map<Species, SpeciesInfoSOA>  InnerSpeciesMapPerEvent;
  typedef std::map<double, InnerSpeciesMapPerEvent> SpeciesMapPerEvent;
  SpeciesMapPerEvent fSpeciesInfoPerEvent;
#endif
  
  std::set<G4double> fTimeToRecord;
  
  int fNEvent; // number of processed events
  double fEdep; // total energy deposition
  G4bool fOutputToRoot; // flag to output to a ROOT file
  G4bool fOutputToXml; // flag to output to a XML file
  G4bool fOutputToCsv; // flag to output to a CSV file (text)

protected:
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

public:
  virtual void Initialize(G4HCofThisEvent*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();
  virtual void DrawAll();
  virtual void PrintAll();
  /** Method used in multithreading mode in order to merge 
      the results*/
  virtual void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* );
  virtual void OutputAndClear();
  
  SpeciesMap GetSpeciesInfo() {return fSpeciesInfoPerTime;}

private:
  G4int fHCID;
  G4THitsMap<G4double>* fEvtMap;
};
#endif
