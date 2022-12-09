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
// J. M. Quesada (July 2009) based on G4GEMCoulombBarrierHE
// Coded strictly according to Furihata's GEM paper 
//
#ifndef G4GEMCoulombBarrier_h
#define G4GEMCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4GEMCoulombBarrier : public G4CoulombBarrier
{
public:
  explicit G4GEMCoulombBarrier(G4int anA, G4int aZ);

  ~G4GEMCoulombBarrier() override = default;

  G4double GetCoulombBarrier(G4int ARes, G4int ZRes, G4double U) const override;

  G4GEMCoulombBarrier(const G4GEMCoulombBarrier & right) = delete;
  const G4GEMCoulombBarrier & operator=
  (const G4GEMCoulombBarrier & right) = delete;

private:
  
  G4double CalcCompoundRadius(G4int ARes) const;

  G4double AejectOneThird;
};
#endif

