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
// $Id: G4CoulombBarrier.cc,v 1.6 2006/06/29 20:28:21 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4CoulombBarrier.hh"
#include "G4HadronicException.hh"
#include <sstream>

G4CoulombBarrier::G4CoulombBarrier(const G4CoulombBarrier & ) : G4VCoulombBarrier()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4CoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4CoulombBarrier & G4CoulombBarrier::operator=(const G4CoulombBarrier & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4CoulombBarrier::operator= meant to not be accessable.");
  return *this;
}

G4bool G4CoulombBarrier::operator==(const G4CoulombBarrier & ) const 
{
  return false;
}

G4bool G4CoulombBarrier::operator!=(const G4CoulombBarrier & ) const 
{
  return true;
}



G4double G4CoulombBarrier::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double U) const 
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
    G4double CompoundRadius = CalcCompoundRadius(static_cast<G4double>(ZRes));
    Barrier = elm_coupling/CompoundRadius * static_cast<G4double>(GetZ())*static_cast<G4double>(ZRes)/
      (std::pow(static_cast<G4double>(GetA()),1./3.) + std::pow(static_cast<G4double>(ARes),1./3.));

    // Barrier penetration coeficient
    G4double K = BarrierPenetrationFactor(ZRes);
		
    Barrier *= K;
		
    Barrier /= (1.0 + std::sqrt(U/(2.0*static_cast<G4double>(ARes))));
  }
  return Barrier;
}



