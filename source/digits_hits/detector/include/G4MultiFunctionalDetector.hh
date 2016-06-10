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
//
// $Id: G4MultiFunctionalDetector.hh 67992 2013-03-13 10:59:57Z gcosmo $
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

