// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReactionKinematics.hh,v 1.2 1999-12-15 14:53:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1995
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ------------ G4ReactionDynamics::TwoBody ------
//                   new inline method BreakupMomentum
//                   new method TwoBodyScattering
//                   by Christian V"olcker (CERN-Munich), August 1997
//                     E-mail: Christian.Volcker@cern.ch
// ************************************************************
//-----------------------------------------------------------------------------

#ifndef G4ReactionKinematics_h
#define G4ReactionKinematics_h 1
 
#include "globals.hh"
#include "Randomize.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"

class G4ReactionKinematics {
public:
   void TwoBodyScattering( const G4DynamicParticle* pIn1, const G4DynamicParticle* pIn2,
               G4DynamicParticle* pOut1, G4DynamicParticle* pOut2);

   inline G4double BreakupMomentum( G4double totalMass, G4double m1, G4double m2);

};

inline G4double G4ReactionKinematics::BreakupMomentum( 
          G4double totalMass, G4double massA, G4double massB){

// is aequivalent to  G4double G4PhaseSpaceDecayChannel::Pmx !!

     G4double m0squared=totalMass*totalMass;
     G4double breakupMomentumSquared=
         (m0squared-(massA+massB)*(massA+massB))*
         (m0squared-(massA-massB)*(massA-massB))/
         (4*m0squared);
    if (breakupMomentumSquared>0) return sqrt(breakupMomentumSquared);
    else       return -1.;
}        


#endif
