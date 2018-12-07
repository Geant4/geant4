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
/// \file runAndEvent/RE02/include/RE02Run.hh
/// \brief Definition of the RE02Run class
//
//
// 

#ifndef RE02Run_h
#define RE02Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include <vector>

//---------------------------------------------------------------------
/// User run class
///
/// (Description) 
///    An example implementation for the multi-functional-detector and 
///   primitive scorers.
///    This RE02Run class has collections which accumulate event information 
///   into run information.
///
/// - constructor
///     gets HitsCollection names, collection IDs and HitsMaps of
///        primitive scorers from the muti-functional detector
///
/// - void RecordEvent(const G4Event*)
///     accumulates HitsMaps over all events into a run HitsMap
///
/// - G4int GetNumberOfHitsMap() const
///     gets the size of the run HitsMap
///
/// - G4THitsMap<G4double>* GetHitsMap(G4int i)
///     gets a run HitsMap of the i-th primitive scorer
///
/// - G4THitsMap<G4double>* GetHitsMap(const G4String& detName, 
///                                  const G4String& colName)
///     gets a run HitsMap with the detName and the colName
///
/// - G4THitsMap<G4double>* GetHitsMap(const G4String& fullName)
///     gets a run HitsMap with the fullName
///
/// - void DumpAllScorer()
///     shows all HitsMap information of this run.
///     This method calls G4THisMap::PrintAll() for individual HitsMap.
//---------------------------------------------------------------------
class RE02Run : public G4Run {

public:
  // constructor and destructor.
  //  vector of multifunctionaldetector name has to given to constructor.
  RE02Run(const std::vector<G4String> mfdName);
  virtual ~RE02Run();

public:
  // virtual method from G4Run. 
  // The method is overriden in this class for scoring.
  virtual void RecordEvent(const G4Event*);
  virtual void Merge(const G4Run*);

  // Access methods for scoring information.
  // - Number of HitsMap for this RUN. 
  //   This is equal to number of collections.
  G4int GetNumberOfHitsMap() const {return fRunMap.size();}
  // - Get HitsMap of this RUN.
  //   by sequential number, by multifucntional name and collection name,
  //   and by collection name with full path.
  G4THitsMap<G4double>* GetHitsMap(G4int i){return fRunMap[i];}
  G4THitsMap<G4double>* GetHitsMap(const G4String& detName, 
                                  const G4String& colName);
  G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);
  // - Dump All HitsMap of this RUN.
  //   This method calls G4THisMap::PrintAll() for individual HitsMap.
  void DumpAllScorer();

private:
  std::vector<G4String> fCollName;
  std::vector<G4int> fCollID;
  std::vector<G4THitsMap<G4double>*> fRunMap;
};

//

#endif
