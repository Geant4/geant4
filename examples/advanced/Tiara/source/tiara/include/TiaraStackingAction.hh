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
//
#ifndef TiaraStackingAction_h
#define TiaraStackingAction_h
#include <map>
#include "globals.hh"
#include "G4UserStackingAction.hh"

class TiaraStackingAction : public G4UserStackingAction
{
public:
  TiaraStackingAction();
  ~TiaraStackingAction();
  virtual G4ClassificationOfNewTrack
  ClassifyNewTrack(const G4Track* aTrack);
  void AddParticleCut(const G4String &particle, G4double cut);
 
private:
  std::map<G4String, G4double> fParticleCut;
  G4double fMinEnergyCut;
};

#endif
