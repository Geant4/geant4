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
// $Id: G4CoulombBarrier.cc,v 1.10 2010-11-15 12:44:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)
//
// 14-11-2007 modified barrier by JMQ (test30) 
// 15-11-2010 V.Ivanchenko use G4Pow and cleanup 

#include "G4CoulombBarrier.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"
#include <sstream>

G4CoulombBarrier::G4CoulombBarrier(): G4VCoulombBarrier(1,0) 
{}

G4CoulombBarrier::G4CoulombBarrier(G4int anA, G4int aZ)
  : G4VCoulombBarrier(anA,aZ) 
{}

G4CoulombBarrier::~G4CoulombBarrier() 
{}

G4double G4CoulombBarrier::BarrierPenetrationFactor(G4double ) const 
{
  return 1.0;
}

G4double G4CoulombBarrier::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double) const 
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
  G4double Barrier = 0.0;
  if (ZRes > ARes || ARes < 1) {
    std::ostringstream errOs;
    errOs << "G4CoulombBarrier::GetCoulombBarrier: ";
    errOs << "Wrong values for ";
    errOs << "residual nucleus A = " << ARes << " ";
    errOs << "and residual nucleus Z = " << ZRes << G4endl;

    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
  if (GetA() == 1 && GetZ() == 0) {
    Barrier = 0.0;   // Neutron Coulomb Barrier is 0
  } else {

    // JMQ: old coulomb barrier commented since it does not agree with Dostrovski's prescription
    // and too low  barriers are obtained (for protons at least)
    // calculation of K penetration factor is correct
    //    G4double CompoundRadius = CalcCompoundRadius(static_cast<G4double>(ZRes));
    //    Barrier = elm_coupling/CompoundRadius * static_cast<G4double>(GetZ())*static_cast<G4double>(ZRes)/
    //      (std::pow(static_cast<G4double>(GetA()),1./3.) + std::pow(static_cast<G4double>(ARes),1./3.));

    ///New coulomb Barrier according to original Dostrovski's paper 
    G4double rho=1.2*fermi; 
    if(GetA()==1 && GetZ()==1){  rho=0.0;}  

    G4double RN=1.5*fermi;  
    // VI cleanup 
    Barrier=elm_coupling*(GetZ()*ZRes)/(RN * G4Pow::GetInstance()->Z13(ARes) + rho);

    // Barrier penetration coeficient
    G4double K = BarrierPenetrationFactor(ZRes);

    Barrier *= K;
		
    // JMQ : the following statement has unknown origin and dimensionally is meaningless( energy divided by mass number in argument of sqrt function). Energy dependence of Coulomb barrier penetrability should be included in proper way (if needed..)
    //   Barrier /= (1.0 + std::sqrt(U/(2.0*static_cast<G4double>(ARes))));
    //
  }
  return Barrier;
}



