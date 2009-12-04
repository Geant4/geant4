//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef ELECTRONRUN_HH
#define ELECTRONRUN_HH

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4THitsMap.hh"
#include <map>

class G4Event;

class ElectronRun : public G4Run {

public:

  ElectronRun(const G4String& detectorName);
  virtual ~ElectronRun();

public:

  virtual void RecordEvent(const G4Event*);
  void DumpData(G4String&) const;

private:

  void Print(const std::vector<G4String>& title,
			 const std::map< G4int, std::vector<G4double> >&out,
			 G4String&) const;  

  std::map<G4int, G4THitsMap<G4double>* > fMap;

};

#endif
