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
// $Id: G4CoulombBarrier.hh,v 1.2 2003/11/03 17:53:06 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
