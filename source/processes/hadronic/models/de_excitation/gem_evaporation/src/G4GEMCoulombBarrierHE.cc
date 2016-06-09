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
// $Id: G4GEMCoulombBarrierHE.cc,v 1.4 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4GEMCoulombBarrierHE.hh"
#include "G4HadronicException.hh"
#include <strstream>

G4GEMCoulombBarrierHE::G4GEMCoulombBarrierHE(const G4GEMCoulombBarrierHE & ) : G4VCoulombBarrier()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GEMCoulombBarrierHE::copy_constructor meant to not be accessable.");
}


const G4GEMCoulombBarrierHE & G4GEMCoulombBarrierHE::operator=(const G4GEMCoulombBarrierHE & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GEMCoulombBarrierHE::operator= meant to not be accessable.");
  return *this;
}

G4bool G4GEMCoulombBarrierHE::operator==(const G4GEMCoulombBarrierHE & ) const 
{
  return false;
}

G4bool G4GEMCoulombBarrierHE::operator!=(const G4GEMCoulombBarrierHE & ) const 
{
  return true;
}



G4double G4GEMCoulombBarrierHE::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double U) const 
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
  G4double Barrier = 0.0;
  if (ZRes > ARes || ARes < 1) {
    char errMessage[1024];
    std::ostrstream errOs(errMessage,1024);
    errOs << "G4GEMCoulombBarrierHE::GetCoulombBarrier: ";
    errOs << "Wrong values for ";
    errOs << "residual nucleus A = " << ARes << " ";
    errOs << "and residual nucleus Z = " << ZRes << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errMessage);
  }
  if (GetZ() == 0) {
    Barrier = 0.0;   // If there is no charge there is neither barrier
  } else {
    G4double CompoundRadius = CalcCompoundRadius(static_cast<G4double>(ARes));
    Barrier = ( elm_coupling * static_cast<G4double>(GetZ()) * static_cast<G4double>(ZRes) )/
      (CompoundRadius+3.75*fermi);
    
    // Barrier penetration coeficient
//    G4double K = BarrierPenetrationFactor(ZRes);
    
//    Barrier *= K;
    
    Barrier /= (1.0 + std::sqrt(U/(2.0*static_cast<G4double>(ARes))));
  }
  return Barrier;
}


G4double G4GEMCoulombBarrierHE::CalcCompoundRadius(const G4double ARes) const
{
    G4double AresOneThird = std::pow(ARes,1.0/3.0);
    G4double AejectOneThird = std::pow(G4double(GetA()),1.0/3.0);

    G4double Result = 1.12*(AresOneThird + AejectOneThird) - 
	0.86*(AresOneThird+AejectOneThird)/(AresOneThird*AejectOneThird);

    return Result*fermi;
}


