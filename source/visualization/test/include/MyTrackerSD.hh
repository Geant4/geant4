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
// $Id: MyTrackerSD.hh,v 1.3 2001-07-11 10:09:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerSD_h
#define MyTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MyTrackerHitsCollection.hh"
#include "MyTrackerHit.hh"
#include "G4Step.hh"

class MyTrackerSD : public G4VSensitiveDetector
{

  public:
      MyTrackerSD(G4String name);
      ~MyTrackerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyTrackerHitsCollection *TrackerCollection;
};




#endif

