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
// $Id: G4ReactionKinematics.hh,v 1.5 2002-12-12 19:18:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1995
//      CERN Geneva Switzerland
//
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
