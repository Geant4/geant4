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
// $Id: G4VPrimitiveSensitivity.hh,v 1.1 2005-09-19 18:40:56 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VPrimitiveSensitivity_h
#define G4VPrimitiveSensitivity_h 1

#include "G4VSensitiveDetector.hh"
#include "G4VHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "G4CollectionNameVector.hh"
#include "G4MultiFunctionalDetector.hh"

// class description:
//
//  This is the base class of the sensitive detector which owns
// only one hits collection and to be registered to G4MultiFunctionalDetector

class G4VPrimitiveSensitivity : public G4VSensitiveDetector
{
  friend class G4MultiFunctionalDetector;

  public: // with description
      G4VPrimitiveSensitivity(G4String);

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*)=0;

      virtual G4int GetIndex(G4Step*);

  public:
      virtual ~G4VPrimitiveSensitivity();

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();

  public:
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      G4MultiFunctionalDetector* detector;

  public:
      virtual G4int GetCollectionID(G4int);
      inline void SetMultiFunctionalDetector(G4MultiFunctionalDetector* d)
      { detector = d; }

};



#endif

