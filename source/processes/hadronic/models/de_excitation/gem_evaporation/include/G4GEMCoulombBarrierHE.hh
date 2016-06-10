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
// $Id: G4GEMCoulombBarrierHE.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4GEMCoulombBarrierHE_h
#define G4GEMCoulombBarrierHE_h 1

#include "G4VCoulombBarrier.hh"
#include "globals.hh"


class G4GEMCoulombBarrierHE : public G4VCoulombBarrier
{
public:

  G4GEMCoulombBarrierHE(G4int anA, G4int aZ);

  ~G4GEMCoulombBarrierHE();

  G4double GetCoulombBarrier(G4int ARes, G4int ZRes, G4double U) const;

private:

  G4GEMCoulombBarrierHE();
  G4GEMCoulombBarrierHE(const G4GEMCoulombBarrierHE & right);
  const G4GEMCoulombBarrierHE & operator=(const G4GEMCoulombBarrierHE & right);
  G4bool operator==(const G4GEMCoulombBarrierHE & right) const;
  G4bool operator!=(const G4GEMCoulombBarrierHE & right) const;
  
public:

  virtual G4double BarrierPenetrationFactor(G4double /*aZ*/) const
  {return 1.0;};

  G4double CalcCompoundRadius(G4int ARes) const;

};
#endif

