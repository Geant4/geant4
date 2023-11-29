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
//

#ifndef G4PSCylinderSurfaceCurrent_h
#define G4PSCylinderSurfaceCurrent_h 1

#include "G4VPrimitivePlotter.hh"
#include "G4THitsMap.hh"

#include "G4Tubs.hh"
#include "G4PSDirectionFlag.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Surface Current.
//  Current version assumes only for G4Tubs shape, and the surface
//  is defined at the inner plane of the tubs.
//   The current is given in the unit of area.
//    e.g.  (Number of tracks)/mm2.
//
// Surface is defined at the  inner surface of the tube.
// Direction                   R    R+dR
//   0  IN || OUT            ->|<-  |      fCurrent_InOut
//   1  IN                   ->|    |      fCurrent_In
//   2  OUT                    |<-  |      fCurrent_Out
//
//
// Created: 2007-03-21  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
// 2020-10-06   Use G4VPrimitivePlotter and fill 1-D histo of kinetic energy (x)
//              vs. surface current * track weight (y)             (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

class G4PSCylinderSurfaceCurrent : public G4VPrimitivePlotter
{
 public: 
  G4PSCylinderSurfaceCurrent(G4String name, G4int direction, G4int depth = 0);
  G4PSCylinderSurfaceCurrent(G4String name, G4int direction,
                             const G4String& unit, G4int depth = 0);
  ~G4PSCylinderSurfaceCurrent() override = default;

  inline void Weighted(G4bool flg = true) { weighted = flg; }
  // Multiply track weight

  inline void DivideByArea(G4bool flg = true) { divideByArea = flg; }
  // Divided by Area.

  void Initialize(G4HCofThisEvent*) override;
  void clear() override;
  void PrintAll() override;

  virtual void SetUnit(const G4String& unit);

 protected: 
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  G4int IsSelectedSurface(G4Step*, G4Tubs*);
  virtual void DefineUnitAndCategory();

 private:
  G4int HCID;
  G4int fDirection;
  G4THitsMap<G4double>* EvtMap;
  G4bool weighted;
  G4bool divideByArea;
};

#endif
