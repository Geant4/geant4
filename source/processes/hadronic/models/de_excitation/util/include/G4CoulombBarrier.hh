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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)
//
// 15-11-2010 V.Ivanchenko cleanup 

#ifndef G4CoulombBarrier_h
#define G4CoulombBarrier_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VCoulombBarrier.hh"
#include "G4HadronicException.hh"

class G4CoulombBarrier : public G4VCoulombBarrier
{

public:

  G4CoulombBarrier();
  G4CoulombBarrier(G4int anA, G4int aZ);
  virtual ~G4CoulombBarrier();

  G4double GetCoulombBarrier(G4int ARes, G4int ZRes, G4double U) const;

private:
  G4CoulombBarrier(const G4CoulombBarrier & right);

  const G4CoulombBarrier & operator=(const G4CoulombBarrier & right);
  G4bool operator==(const G4CoulombBarrier & right) const;
  G4bool operator!=(const G4CoulombBarrier & right) const;
  
  virtual G4double BarrierPenetrationFactor(G4double ) const;

  inline G4double CalcCompoundRadius(const G4double ZRes) const 
  {
    return 2.173*CLHEP::fermi*(1.0+0.006103*static_cast<G4double>(GetZ())*ZRes)/
      (1.0+0.009443*static_cast<G4double>(GetZ())*ZRes);
  }
};
#endif
