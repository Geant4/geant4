// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReactionKinematics.cc,v 1.2 1999-12-15 14:53:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      ------------ G4ReactionDynamics::TwoBody ------
//                     new method TwoBodyScattering
//                   by Christian V"olcker (CERN-Munich), September 1997
//                     E-mail: Christian.Volcker@cern.ch
// ************************************************************
//-----------------------------------------------------------------------------

#include "G4ReactionKinematics.hh"

// ************************************************************
void G4ReactionKinematics::TwoBodyScattering(
           const G4DynamicParticle* pIn1, const G4DynamicParticle* pIn2,
           G4DynamicParticle* pOut1, G4DynamicParticle* pOut2)
// ************************************************************
{           
// initial particles:

// - total invariant mass
   G4LorentzVector sumIn(pIn1->Get4Momentum()+pIn2->Get4Momentum());
   G4double invariantMass=sumIn.mag();

// - beta of center-of-mass system
   G4ThreeVector betaCMS=sumIn.boostVector();

// final particles:

// - get final particle masses
   G4double massOut1=pOut1->GetMass();
   G4double massOut2=pOut2->GetMass();

// - calculate breakup momentum:
   G4double breakupMomentum=BreakupMomentum(invariantMass, massOut1, massOut2);

// - random decay angle
   G4double theta=RandFlat::shoot(HepDouble(0.),HepDouble(pi));  // isotropic decay angle theta
   G4double phi  =RandFlat::shoot(HepDouble(0.),HepDouble(twopi));  // isotropic decay angle phi

// - setup LorentzVectors
   G4double pz=cos(theta)*breakupMomentum;
   G4double px=sin(theta)*cos(phi)*breakupMomentum;
   G4double py=sin(theta)*sin(phi)*breakupMomentum;
   
   G4double breakupMomentumSquared=breakupMomentum*breakupMomentum;
   G4double energy1=sqrt(breakupMomentumSquared+massOut1*massOut1);
   G4double energy2=sqrt(breakupMomentumSquared+massOut2*massOut2);

   G4LorentzVector lorentz1(px, py, pz, energy1);
   G4LorentzVector lorentz2(px, py, pz, energy2);

// - back into lab system

   lorentz1.boost(betaCMS);
   lorentz2.boost(betaCMS);

// fill in new particles:

   pOut1->Set4Momentum(lorentz1);
   pOut2->Set4Momentum(lorentz2);

   return;
}
