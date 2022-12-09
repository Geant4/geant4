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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)
//
// 14-11-2007 modified barrier by JMQ (test30) 
// 15-11-2010 V.Ivanchenko use G4Pow and cleanup 

#include "G4CoulombBarrier.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearRadii.hh"

G4CoulombBarrier::G4CoulombBarrier(G4int A, G4int Z)
  : G4VCoulombBarrier(A, Z)
{
  factor = CLHEP::elm_coupling*Z;
  SetParameters(0.4*G4NuclearRadii::RadiusCB(Z, A), 1.5*CLHEP::fermi);
}

G4double G4CoulombBarrier::GetCoulombBarrier(
         G4int ARes, G4int ZRes, G4double U) const 
{
  if(0 == theZ) { return 0.0; }
  G4double cb = factor*ZRes/(G4NuclearRadii::RadiusCB(ZRes,ARes) + theRho);
  if(U > 0.0) { cb /= (1.0 + std::sqrt( U/((2*ARes)*CLHEP::MeV) )); }
  return cb;
}

G4double G4CoulombBarrier::BarrierPenetrationFactor(G4int aZ) const 
{
  // Data comes from 
  // Dostrovsky, Fraenkel and Friedlander
  // Physical Review, vol 116, num. 3 1959
  // 
  // const G4int size = 5;
  // const G4double Zlist[size] = {10.0, 20.0, 30.0, 50.0, 70.0};
  // const G4double Kprot[size] = {0.42, 0.58, 0.68, 0.77, 0.80};
  // 
  G4double res = 1.0;
  if(theZ == 1) {
    res = (aZ >= 70) ? 0.80 :
    (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
    res += 0.06*(theA - 1);

  } else if(theZ == 2 && theA <= 4) {
    res = (aZ >= 70) ? 0.98 :
    (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
    res += 0.12*(4 - theA);
  }
  return res;
}
