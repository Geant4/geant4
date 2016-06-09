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
// $Id: G4PSPassageTrackLength.hh,v 1.2 2005/11/19 00:44:00 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4PSPassageTrackLength_h
#define G4PSPassageTrackLength_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only track length.
//   The tracks which passed a geometry is taken into account.
// 
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

class G4PSPassageTrackLength : public G4VPrimitiveScorer
{
 
  public: // with description
      G4PSPassageTrackLength(G4String name, G4int depth=0);
      virtual ~G4PSPassageTrackLength();

      inline void Weighted(G4bool flg=true) { weighted = flg; }
      // Multiply track weight

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      G4bool IsPassed(G4Step*);

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      G4int HCID;
      G4int fCurrentTrkID;
      G4double fTrackLength;
      G4THitsMap<G4double>* EvtMap;
      G4bool weighted;

};

#endif

