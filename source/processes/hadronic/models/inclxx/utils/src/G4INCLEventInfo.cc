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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLEventInfo.cc
 * \brief Simple container for output of event results.
 *
 * Contains the results of an INCL cascade.
 *
 * \date 21 January 2011
 * \author Davide Mancusi
 */

#include "G4INCLEventInfo.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLParticleTable.hh"
#include <cmath>

namespace G4INCL {

  Int_t EventInfo::eventNumber = 0;

#ifdef INCL_INVERSE_KINEMATICS
  void EventInfo::fillInverseKinematics(const Double_t gamma) {
    const Double_t beta = std::sqrt(1.-1./(gamma*gamma));
    for(Int_t i=0; i<nParticles; ++i) {
      Double_t mass;
      if(A[i]>0) {
        mass = ParticleTable::getTableMass(A[i],Z[i]);
      } else if(origin[i]==-1) { // cascade particles with A=0, must be pions
        if(Z[i]==1)
          mass = ParticleTable::getTableParticleMass(PiPlus);
        else if(Z[i]==0)
          mass = ParticleTable::getTableParticleMass(PiZero);
        else
          mass = ParticleTable::getTableParticleMass(PiMinus);
      } else // gamma rays
        mass = 0.;

      const Double_t ETot = EKin[i] + mass;
      const Double_t ETotPrime = gamma*(ETot - beta*pz[i]);
      /* Using the invariant mass here avoids negative kinetic energies with
       * particles produced by the de-excitation models, which do not
       * necessarily use the same mass look-up tables as INCL.
       */
      Double_t invariantMass;
      if(A[i]>0 || origin[i]==-1) { // massive particles
        invariantMass = std::sqrt(ETot*ETot - px[i]*px[i] - py[i]*py[i] - pz[i]*pz[i]);
      } else { // gamma rays
        invariantMass = 0.;
      }
      EKinPrime[i] = ETotPrime - invariantMass;
      pzPrime[i] = -gamma*(pz[i] - beta*ETot);
      const Double_t pPrime = std::sqrt(px[i]*px[i] + py[i]*py[i] + pzPrime[i]*pzPrime[i]);
      const Double_t cosThetaPrime = pzPrime[i]/pPrime;
      if(cosThetaPrime>=1.)
        thetaPrime[i] = 0.;
      else if(cosThetaPrime<=-1.)
        thetaPrime[i] = 180.;
      else
        thetaPrime[i] = 180.*std::acos(cosThetaPrime)/Math::pi;
    }
  }
#endif // INCL_INVERSE_KINEMATICS
}

