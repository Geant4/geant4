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

#ifndef G4HumanPhantomSDBreast_h
#define G4HumanPhantomSDBreast_h 1

#include "G4VSensitiveDetector.hh"


class G4Step;
class G4HCofThisEvent;

class G4HumanPhantomSDBreast : public G4VSensitiveDetector
{
  public:
      G4HumanPhantomSDBreast(G4String);
     ~G4HumanPhantomSDBreast();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  private:
  // G4HumanPhantomHitsCollection* trackerBreastCollection;

};


#endif

