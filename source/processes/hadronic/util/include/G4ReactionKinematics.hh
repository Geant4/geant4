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
    if (breakupMomentumSquared>0) return std::sqrt(breakupMomentumSquared);
    else       return -1.;
}        


#endif
