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

// ====================================================================
//
//   XXMuonSD.hh
//   $Id: XXMuonSD.hh,v 1.1 2002-04-29 20:44:41 asaim Exp $
//
// ====================================================================
#ifndef XX_MUON_SD_H
#define XX_MUON_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "XXMuonHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class XXMuonSD : public G4VSensitiveDetector {
private:
  XXMuonHitsCollection* hitCollection;

public:
  XXMuonSD(G4String name);
  ~XXMuonSD();
  
  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory* ROhist);
  virtual void EndOfEvent(G4HCofThisEvent* HCE);
  virtual void clear();
  virtual void DrawAll();
  virtual void PrintAll();
  
};

#endif
