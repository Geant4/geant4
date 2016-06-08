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
#ifndef G4Absorber_hh
#define G4Absorber_hh

#include "globals.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;

class G4Absorber
{
public:
  G4Absorber(G4double cutOnP);
  ~G4Absorber();

  G4bool WillBeAbsorbed(const G4KineticTrack & kt);
  G4bool Absorb(G4KineticTrack & kt, G4KineticTrackVector & tgt);
  G4KineticTrackVector * GetAbsorbers();
  G4KineticTrackVector * GetProducts();
  G4bool FindAbsorbers(G4KineticTrack & kt, G4KineticTrackVector & tgt);
  G4bool FindProducts(G4KineticTrack & kt);


private:
  G4double theCutOnP;
  G4KineticTrackVector * theAbsorbers;
  G4KineticTrackVector * theProducts;

  G4ThreeVector GetRandomDirection();


};


inline G4KineticTrackVector * G4Absorber::GetAbsorbers()
{
  return theAbsorbers;
}

inline G4KineticTrackVector * G4Absorber::GetProducts()
{
  return theProducts;
}

#endif


