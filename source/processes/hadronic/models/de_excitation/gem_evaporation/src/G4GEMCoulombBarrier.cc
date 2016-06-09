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
// J. M. Quesada (July 2009):  New class based on G4GEMCoulombBarrierHE
// Coded strictly according to Furihata's GEM paper 
// NEW:effective decrease  of barrier with E* (Barashenkov) has been added
//
#include "G4GEMCoulombBarrier.hh"
#include "G4HadronicException.hh"
#include <sstream>

G4GEMCoulombBarrier::G4GEMCoulombBarrier(const G4GEMCoulombBarrier & ) : G4VCoulombBarrier()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GEMCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4GEMCoulombBarrier & G4GEMCoulombBarrier::operator=(const G4GEMCoulombBarrier & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GEMCoulombBarrier::operator= meant to not be accessable.");
  return *this;
}

G4bool G4GEMCoulombBarrier::operator==(const G4GEMCoulombBarrier & ) const 
{
  return false;
}

G4bool G4GEMCoulombBarrier::operator!=(const G4GEMCoulombBarrier & ) const 
{
  return true;
}



G4double G4GEMCoulombBarrier::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double U) const 
// Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
  G4double Barrier = 0.0;
  if (ZRes > ARes || ARes < 1) {
    std::ostringstream errOs;
    errOs << "G4GEMCoulombBarrier::GetCoulombBarrier: ";
    errOs << "Wrong values for ";
    errOs << "residual nucleus A = " << ARes << " ";
    errOs << "and residual nucleus Z = " << ZRes << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
  if (GetZ() == 0) {
    Barrier = 0.0;   // If there is no charge there is neither barrier
  } else 
    {
      G4double CompoundRadius = CalcCompoundRadius(static_cast<G4double>(ARes));
      Barrier = ( elm_coupling * static_cast<G4double>(GetZ()) * static_cast<G4double>(ZRes) )/CompoundRadius;
      
      // Barrier penetration coeficient
      G4double K=1.;
      if(GetA() <= 4) K = BarrierPenetrationFactor(ZRes);
      
      Barrier *= K;
      
    }
  //JMQ 200709 effective decrease  of barrier with E* (Barashenkov)
  // (not inclued in original Furihata's formulation)
  Barrier /= (1.0 + std::sqrt(U/(2.0*static_cast<G4double>(ARes))));
  // JMQ end test
  return Barrier;
}


G4double G4GEMCoulombBarrier::CalcCompoundRadius(const G4double ARes) const
{      
    G4double AresOneThird = std::pow(ARes,1.0/3.0);
    G4double AejectOneThird = std::pow(G4double(GetA()),1.0/3.0);

    if(GetA()==1){
    G4double Rd=1.7* AresOneThird;
    return Rd*fermi;
   } else if (GetA()==2 || GetA()==3 || GetA()==4){
    G4double Rd=1.7* AresOneThird;
    G4double Rj=1.2;
    return (Rd+Rj)*fermi;
   } else {
    G4double Result = 1.12*(AresOneThird + AejectOneThird) - 
      0.86*(AresOneThird+AejectOneThird)/(AresOneThird*AejectOneThird)+3.75;
    return Result*fermi;}
}


