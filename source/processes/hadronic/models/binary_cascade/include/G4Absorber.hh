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
  // hide copy ctor, =, == and != operators
  G4Absorber(const  G4Absorber &right);
  const G4Absorber & operator=(const G4Absorber & right);
  int operator==(const G4Absorber & right) const;
  int operator!=(const G4Absorber & right) const;


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


