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
// $Id: G4CoulombBarrier.cc 100690 2016-10-31 11:25:43Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)
//
// 14-11-2007 modified barrier by JMQ (test30) 
// 15-11-2010 V.Ivanchenko use G4Pow and cleanup 

#include "G4CoulombBarrier.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

G4CoulombBarrier::G4CoulombBarrier(G4int A, G4int Z)
  : G4VCoulombBarrier(A, Z) 
{
  g4calc = G4Pow::GetInstance();
  if(Z > 0) {
    G4double rho = 1.2*CLHEP::fermi; 
    G4double r0  = 1.5*CLHEP::fermi; 
    if(1 == A) {
      rho = 0.0;
    } else if(A <= 3) {
      rho = 0.8*CLHEP::fermi; 
      r0  = 1.7*CLHEP::fermi;
    } else {
      r0  = 1.7*CLHEP::fermi;
    }
    SetParameters(rho, r0);
  }
}

G4CoulombBarrier::~G4CoulombBarrier() 
{}

G4double G4CoulombBarrier::GetCoulombBarrier(G4int ARes, G4int ZRes, G4double) const 
{
  return CLHEP::elm_coupling*(GetZ()*ZRes)/(GetR0()*g4calc->Z13(ARes) + GetRho());
}

G4double G4CoulombBarrier::BarrierPenetrationFactor(G4int aZ) const 
{
  G4double res = 1.0;
  if(GetZ() == 1) {
    res = (aZ >= 70) ? 0.80 :
    (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
    res += 0.06*(GetA() - 1);

  } else if(GetZ() == 2 && GetA() <= 4) {
    res = (aZ >= 70) ? 0.98 :
    (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
    res += 0.12*(4 - GetA());
  }
  return res;
}
