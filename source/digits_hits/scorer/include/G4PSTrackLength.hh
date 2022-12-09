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

#ifndef G4PSTrackLength_h
#define G4PSTrackLength_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring track length.
//
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// Modified: 2007-02-02 Tsukasa ASO, Add MultiplyKineticEnergy()
//                                  and DivideByVelocity().
//           2010-07-22   Introduce Unit specification.
//
//
///////////////////////////////////////////////////////////////////////////////

class G4PSTrackLength : public G4VPrimitiveScorer
{
 public:
  G4PSTrackLength(G4String name, G4int depth = 0);
  G4PSTrackLength(G4String name, const G4String& unit, G4int depth = 0);
  ~G4PSTrackLength() override = default;

  inline void Weighted(G4bool flg = true) { weighted = flg; }
  // Multiply track weight

  void MultiplyKineticEnergy(G4bool flg = true);
  // Multiply Kinetic Energy

  void DivideByVelocity(G4bool flg = true);
  // Divide by velocity

  void Initialize(G4HCofThisEvent*) override;
  void clear() override;
  void PrintAll() override;

  virtual void SetUnit(const G4String& unit);

 protected:
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  virtual void DefineUnitAndCategory();

 private:
  G4int HCID;
  G4THitsMap<G4double>* EvtMap;
  G4bool weighted;
  G4bool multiplyKinE;
  G4bool divideByVelocity;
};
#endif
