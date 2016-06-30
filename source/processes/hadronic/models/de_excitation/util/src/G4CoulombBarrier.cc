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
// $Id: G4CoulombBarrier.cc 97097 2016-05-25 07:50:31Z gcosmo $
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

G4CoulombBarrier::G4CoulombBarrier(G4int anA, G4int aZ)
  : G4VCoulombBarrier(anA,aZ) 
{
  g4pow = G4Pow::GetInstance();
}

G4CoulombBarrier::~G4CoulombBarrier() 
{}

G4double G4CoulombBarrier::GetCoulombBarrier(G4int ARes, G4int ZRes, G4double) const 
{
  G4double Barrier = 0.0;
  if (GetZ() > 0 && ZRes > 0) {
    // JMQ: old coulomb barrier commented since it does not agree with 
    //      Dostrovski's prescription
    //      and too low  barriers are obtained (for protons at least)

    ///New coulomb Barrier according to original Dostrovski's paper 
    static const G4double rho0 = 1.2*CLHEP::fermi; 
    G4double rho = (GetA()==1 && GetZ()==1) ?  0.0 : rho0;  

    static const G4double RN = 1.5*fermi;  
    // VI cleanup 
    Barrier = CLHEP::elm_coupling*(GetZ()*ZRes)/(RN * g4pow->Z13(ARes) + rho);

    // Barrier penetration coeficient
    Barrier *= BarrierPenetrationFactor(ZRes);

  }
  return Barrier;
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
