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
// $Id: G4PSTermination.hh,v 1.1 2007-08-14 21:23:51 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PSTermination_h
#define G4PSTermination_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"


//////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Number of Steps in the cell.
// 
// Created: 2007-2-2  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////


class G4PSTermination : public G4VPrimitiveScorer
{
 
 public: // with description
      G4PSTermination(G4String name, G4int depth=0);

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

  public:
      virtual ~G4PSTermination();

      inline void Weighted(G4bool flg=true) { weighted = flg; }
      // Multiply track weight

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();

  public:
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      G4int HCID;
      G4THitsMap<G4double>* EvtMap;
      G4bool weighted;

  public:


};



#endif
