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
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
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
#include "G4INCLParticle.hh"
#include <cmath>

namespace G4INCL {

  G4ThreadLocal Int_t EventInfo::eventNumber = 0;

  void EventInfo::fillInverseKinematics(const Double_t gamma) {
    const Double_t beta = std::sqrt(1.-1./(gamma*gamma));
    for(Int_t i=0; i<nParticles; ++i) {
      // determine the particle mass from the kinetic energy and the momentum;
      // this ensures consistency with the masses uses by the models
      Double_t mass;
      if(EKin[i]>0.) {
        mass = std::max(
                        0.5 * (px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]-EKin[i]*EKin[i]) / EKin[i],
                        0.0);
      } else {
        INCL_WARN("Particle with null kinetic energy in fillInverseKinematics, cannot determine its mass:\n"
                  << "  A=" << A[i] << ", Z=" << Z[i] << ", S=" << S[i] << '\n'
                  << "  EKin=" << EKin[i] << ", px=" << px[i] << ", py=" << py[i] << ", pz=" << pz[i] << '\n'
                  << "  Falling back to the mass from the INCL ParticleTable" << '\n');
        mass = ParticleTable::getRealMass(A[i], Z[i], S[i]);
      }

      const Double_t ETot = EKin[i] + mass;
      const Double_t ETotPrime = gamma*(ETot - beta*pz[i]);
      EKinPrime[i] = ETotPrime - mass;
      pzPrime[i] = -gamma*(pz[i] - beta*ETot);
      const Double_t pPrime = std::sqrt(px[i]*px[i] + py[i]*py[i] + pzPrime[i]*pzPrime[i]);
      const Double_t cosThetaPrime = (pPrime>0.) ? (pzPrime[i]/pPrime) : 1.;
      if(cosThetaPrime>=1.)
        thetaPrime[i] = 0.;
      else if(cosThetaPrime<=-1.)
        thetaPrime[i] = 180.;
      else
        thetaPrime[i] = Math::toDegrees(Math::arcCos(cosThetaPrime));
    }
  }

  void EventInfo::remnantToParticle(const G4int remnantIndex) {
    
    INCL_DEBUG("remnantToParticle function used\n");
    
    A[nParticles] = ARem[remnantIndex];
    Z[nParticles] = ZRem[remnantIndex];
    S[nParticles] = SRem[remnantIndex];
    
	ParticleSpecies pt(A[nParticles],Z[nParticles],S[nParticles]);
	PDGCode[nParticles] = pt.getPDGCode();
	
    ParticleBias[nParticles] = Particle::getTotalBias();
    emissionTime[nParticles] = stoppingTime;

    px[nParticles] = pxRem[remnantIndex];
    py[nParticles] = pyRem[remnantIndex];
    pz[nParticles] = pzRem[remnantIndex];

    const G4double plab = std::sqrt(pxRem[remnantIndex]*pxRem[remnantIndex]
                                  +pyRem[remnantIndex]*pyRem[remnantIndex]
                                  +pzRem[remnantIndex]*pzRem[remnantIndex]);
    G4double pznorm = pzRem[remnantIndex]/plab;
    if(pznorm>1.)
      pznorm = 1.;
    else if(pznorm<-1.)
      pznorm = -1.;
    theta[nParticles] = Math::toDegrees(Math::arcCos(pznorm));
    phi[nParticles] = Math::toDegrees(std::atan2(pyRem[remnantIndex],pxRem[remnantIndex]));

    EKin[nParticles] = EKinRem[remnantIndex];
    origin[nParticles] = -1; // Origin: cascade
    parentResonancePDGCode[nParticles] = 0;  // No parent resonance
    parentResonanceID[nParticles] = 0;       // No parent resonance
    history.push_back(""); // history
    nParticles++;
// assert(history.size()==(unsigned int)nParticles);
  }
}

