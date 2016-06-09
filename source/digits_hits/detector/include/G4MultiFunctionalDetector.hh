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
// $Id: G4MultiFunctionalDetector.hh,v 1.2 2005/11/16 22:59:01 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4MultiFunctionalDetector_h
#define G4MultiFunctionalDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4VHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "G4CollectionNameVector.hh"
#include <vector>
class G4VPrimitiveScorer;

// class description:
//
//  This is the base class of the sensitive detector which owns
// one or more G4VPrimitiveScorer class objects.

class G4MultiFunctionalDetector : public G4VSensitiveDetector
{

  public: // with description
      G4MultiFunctionalDetector(G4String);

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

      std::vector<G4VPrimitiveScorer*> primitives;

  public: // with description
      G4bool RegisterPrimitive(G4VPrimitiveScorer*);
      G4bool RemovePrimitive(G4VPrimitiveScorer*);
      inline G4int GetNumberOfPrimitives() const
      { return primitives.size(); }
      G4VPrimitiveScorer* GetPrimitive(G4int id) const
      { return primitives[id]; }

  public:
      virtual ~G4MultiFunctionalDetector();

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();

  public:
      virtual void DrawAll();
      virtual void PrintAll();


};



#endif

