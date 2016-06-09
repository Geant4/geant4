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
// $Id: G4GEMCoulombBarrierHE.cc,v 1.2 2002/12/12 19:17:20 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4GEMCoulombBarrierHE.hh"
#include "g4std/strstream"

G4GEMCoulombBarrierHE::G4GEMCoulombBarrierHE(const G4GEMCoulombBarrierHE & right)
{
  G4Exception("G4GEMCoulombBarrierHE::copy_constructor meant to not be accessable.");
}


const G4GEMCoulombBarrierHE & G4GEMCoulombBarrierHE::operator=(const G4GEMCoulombBarrierHE & right)
{
  G4Exception("G4GEMCoulombBarrierHE::operator= meant to not be accessable.");
  return *this;
}

G4bool G4GEMCoulombBarrierHE::operator==(const G4GEMCoulombBarrierHE & right) const 
{
  return false;
}

G4bool G4GEMCoulombBarrierHE::operator!=(const G4GEMCoulombBarrierHE & right) const 
{
  return true;
}



G4double G4GEMCoulombBarrierHE::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double U) const 
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
  G4double Barrier = 0.0;
  if (ZRes > ARes || ARes < 1) {
    char errMessage[1024];
    G4std::ostrstream errOs(errMessage,1024);
    errOs << "G4GEMCoulombBarrierHE::GetCoulombBarrier: ";
    errOs << "Wrong values for ";
    errOs << "residual nucleus A = " << ARes << " ";
    errOs << "and residual nucleus Z = " << ZRes << G4endl;
    G4Exception(errMessage);
  }
  if (GetZ() == 0) {
    Barrier = 0.0;   // If there is no charge there is neither barrier
  } else {
    G4double CompoundRadius = CalcCompoundRadius(G4double(ARes));
    Barrier = ( elm_coupling *G4double(GetZ()) * G4double(ZRes) )/
      (CompoundRadius+3.75*fermi);
    
    // Barrier penetration coeficient
//    G4double K = BarrierPenetrationFactor(ZRes);
    
//    Barrier *= K;
    
    Barrier /= (1.0 + sqrt(U/(2.0*G4double(ARes))));
  }
  return Barrier;
}


G4double G4GEMCoulombBarrierHE::CalcCompoundRadius(const G4double ARes) const
{
    G4double AresOneThird = pow(ARes,1.0/3.0);
    G4double AejectOneThird = pow(this->GetA(),1.0/3.0);

    G4double Result = 1.12*(AresOneThird + AejectOneThird) - 
	0.86*(AresOneThird+AejectOneThird)/(AresOneThird*AejectOneThird);

    return Result*fermi;
}


