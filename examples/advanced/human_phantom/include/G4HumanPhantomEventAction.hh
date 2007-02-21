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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
// 
#ifndef G4HumanPhantomEventAction_h
#define G4HumanPhantomEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>

class G4Event;
class G4HumanPhantomEventAction : public G4UserEventAction
{
public:
  G4HumanPhantomEventAction();
  ~G4HumanPhantomEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private:
  void Fill(G4String bodypartName, G4double energyDeposit);
  void totalEventEnergyDeposit();
 
  G4int hitCollectionID; 
  std::map<std::string,G4double> energyTotal;
  G4String bodypartName;
};
#endif

    
