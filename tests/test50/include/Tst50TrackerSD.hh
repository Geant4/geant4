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
// $Id: Tst50TrackerSD.hh,v 1.4 2003-07-03 13:43:10 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG      first implemntation
// -------------------------------------------------------------------

#ifndef Tst50TrackerSD_h
#define Tst50TrackerSD_h 1

#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "Tst50TrackerHit.hh"

class G4Step;
class G4HCofThisEvent;

class Tst50DetectorConstruction;

class Tst50TrackerSD : public G4VSensitiveDetector
{
public:
  Tst50TrackerSD(G4String, Tst50DetectorConstruction*);
  ~Tst50TrackerSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  Tst50TrackerHitsCollection* tst50Collection;
  G4int*                   hitID;
  Tst50DetectorConstruction* Detector;
};
#endif

