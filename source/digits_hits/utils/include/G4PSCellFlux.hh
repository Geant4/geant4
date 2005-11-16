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
// $Id: G4PSCellFlux.hh,v 1.1 2005-11-16 23:12:42 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PSCellFlux_h
#define G4PSCellFlux_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only cell flux.
//   The Cell Flux is defined by  a track length divided by a geometry
//   volume, where all of the tracks in the geometry are taken 
//  into account. 
//    Please use G4PSPassageCellFlux, if you want to score only 
//  tracks which passed through the volume.
//
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////


class G4PSCellFlux : public G4VPrimitiveScorer
{
   public: // with description
      G4PSCellFlux(G4String name, G4int depth=0);
      virtual ~G4PSCellFlux();

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      G4int HCID;
      G4THitsMap<G4double>* EvtMap;
};
#endif

