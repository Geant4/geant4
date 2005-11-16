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
// $Id: G4PSSphereSurfaceCurrent.hh,v 1.1 2005-11-16 23:12:42 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PSSphereSurfaceCurrent_h
#define G4PSSphereSurfaceCurrent_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

#include "G4Sphere.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Surface Current.
//  Current version assumes only for G4Sphere shape. 
//
// Surface is defined  at the inside of sphere.
// Direction                  -Rmin   +Rmax
//   0  IN || OUT            ->|<-     |
//   1  IN                   ->|       |
//   2  OUT                    |<-     |
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

enum { fCurrent_InOut, fCurrent_In, fCurrent_Out };

class G4PSSphereSurfaceCurrent : public G4VPrimitiveScorer
{
 
  public: // with description
      G4PSSphereSurfaceCurrent(G4String name, G4int direction, G4int depth=0);
      virtual ~G4PSSphereSurfaceCurrent();

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      G4int IsSelectedSurface(G4Step*,G4Sphere*);

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      G4int  HCID;
      G4int  fDirection;
      G4THitsMap<G4double>* EvtMap;

};

#endif

