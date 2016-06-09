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
// $Id: G4CoulombBarrier.hh,v 1.4 2006/06/29 20:27:39 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4CoulombBarrier_h
#define G4CoulombBarrier_h 1

#include "G4VCoulombBarrier.hh"
#include "G4HadronicException.hh"
#include "globals.hh"


class G4CoulombBarrier : public G4VCoulombBarrier
{
public:
  G4CoulombBarrier() : G4VCoulombBarrier(1,0) {};
  G4CoulombBarrier(const G4int anA,const G4int aZ) :
    G4VCoulombBarrier(anA,aZ) {};
  ~G4CoulombBarrier() {};

private:
  G4CoulombBarrier(const G4CoulombBarrier & right);

  const G4CoulombBarrier & operator=(const G4CoulombBarrier & right);
  G4bool operator==(const G4CoulombBarrier & right) const;
  G4bool operator!=(const G4CoulombBarrier & right) const;
  
public:
  G4double GetCoulombBarrier(const G4int ARes, const G4int ZRes, 
			     const G4double U) const;


private:

  virtual G4double BarrierPenetrationFactor(const G4double ) const {return 1.0;}

  virtual G4double CalcCompoundRadius(const G4double ZRes) const 
  {
    return 2.173*fermi*(1.0+0.006103*static_cast<G4double>(GetZ())*ZRes)/
      (1.0+0.009443*static_cast<G4double>(GetZ())*ZRes);
  }
};
#endif
