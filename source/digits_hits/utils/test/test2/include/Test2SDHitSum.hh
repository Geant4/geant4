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

#ifndef Test2SDHitSum_h
#define Test2SDHitSum_h 1

#include "Test2PhantomHit.hh"
#include "G4THitsMap.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class Test2SDHitSum {

public:
  Test2SDHitSum(const G4String& detName,
		std::vector<G4String>& hcnameVec);
  ~Test2SDHitSum();

  void Analyze(Test2PhantomHit* hc);
  G4double GetTotal(G4int i) const;

  void     DumpQuantitiesToFile();
  void     DumpQuantitiyToFile(G4int i);

private:
  G4double GetTotal(const G4THitsMap<G4double>* map) const;
  void     DumpQuantitiyToFile(const G4THitsMap<G4double> *map, 
			       G4String& fileName);

private:
  G4String              DTName;
  std::vector<G4String> HCName;
  std::vector<G4THitsMap<G4double>*> MapSum;
  

};
#endif

