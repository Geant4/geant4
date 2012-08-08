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
// INCL++ revision: v5.1.2
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLEventInfo.hh
 * \brief Simple container for output of event results.
 *
 * Contains the results of an INCL cascade.
 *
 * Created on: 21 January 2011
 *     Author: Davide Mancusi
 */

#ifndef G4INCLEVENTINFO_HH
#define G4INCLEVENTINFO_HH 1

#include "G4INCLParticleType.hh"
#ifdef INCL_ROOT_USE
#include <Rtypes.h>
#endif
#include <string>
#include <vector>
#include <algorithm>

namespace G4INCL {
#ifndef INCL_ROOT_USE
    typedef G4int Int_t;
    typedef short Short_t;
    typedef G4float Float_t;
    typedef G4bool Bool_t;
#endif

    struct EventInfo {
      EventInfo() :
        projectileType(UnknownParticle),
        At(0), Zt(0), Ap(0), Zp(0),
        Ep(0.),
        impactParameter(0.0), nCollisions(0), stoppingTime(0.0),
        EBalance(0.0), pLongBalance(0.0), pTransBalance(0.0),
        nCascadeParticles(0), nRemnants(0), nParticles(0),
        transparent(true),
        forcedCompoundNucleus(false),
        nucleonAbsorption(false), pionAbsorption(false), nDecays(0),
        nBlockedCollisions(0), nBlockedDecays(0),
        effectiveImpactParameter(0.0),
        deltasInside(false),
        forcedDeltasInside(false),
        forcedDeltasOutside(false),
        clusterDecay(false),
        firstCollisionTime(0.),
        firstCollisionXSec(0.),
        nReflectionAvatars(0),
        nCollisionAvatars(0),
        nDecayAvatars(0),
        nUnmergedSpectators(0)
      {
        std::fill_n(ARem, maxSizeRemnants, 0);
        std::fill_n(ZRem, maxSizeRemnants, 0);
        std::fill_n(EStarRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(JRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(EKinRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(pxRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(pyRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(pzRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(thetaRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(phiRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(jxRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(jyRem, maxSizeRemnants, ((Float_t)0.));
        std::fill_n(jzRem, maxSizeRemnants, ((Float_t)0.));

        std::fill_n(A, maxSizeParticles, 0);
        std::fill_n(Z, maxSizeParticles, 0);
        std::fill_n(emissionTime, maxSizeParticles, ((Float_t)0.));
        std::fill_n(EKin, maxSizeParticles, ((Float_t)0.));
        std::fill_n(px, maxSizeParticles, ((Float_t)0.));
        std::fill_n(py, maxSizeParticles, ((Float_t)0.));
        std::fill_n(pz, maxSizeParticles, ((Float_t)0.));
        std::fill_n(theta, maxSizeParticles, ((Float_t)0.));
        std::fill_n(phi, maxSizeParticles, ((Float_t)0.));
        std::fill_n(origin, maxSizeParticles, 0);
      };

      /** \brief Number of the event */
      static Int_t eventNumber;
      /** \brief Protjectile particle type */
      ParticleType projectileType;

      /** \brief Mass number of the target nucleus */
      Short_t At;
      /** \brief Charge number of the target nucleus */
      Short_t Zt;

      /** \brief Mass number of the projectile nucleus */
      Short_t Ap;
      /** \brief Charge number of the projectile nucleus */
      Short_t Zp;
      /** \brief Projectile kinetic energy given as input */
      Float_t Ep;

      /** \brief Impact parameter [fm] */
      Float_t impactParameter;
      /** \brief Number of accepted two-body collisions */
      Int_t nCollisions;
      /** \brief Cascade stopping time [fm/c] */
      Float_t stoppingTime;

      /** \brief Energy-conservation balance [MeV] */
      Float_t EBalance;
      /** \brief Longitudinal momentum-conservation balance [MeV/c] */
      Float_t pLongBalance;
      /** \brief Transverse momentum-conservation balance [MeV/c] */
      Float_t pTransBalance;

      /** \brief Number of cascade particles */
      Short_t nCascadeParticles;
      /** \brief Number of remnants */
      Int_t nRemnants;
      /** \brief Total number of emitted particles */
      Int_t nParticles;

      /** \brief True if the event is transparent */
      Bool_t transparent;
      /** \brief True if the event is a forced CN */
      Bool_t forcedCompoundNucleus;
      /** \brief True if the event is absorption */
      Bool_t nucleonAbsorption;
      /** \brief True if the event is absorption */
      Bool_t pionAbsorption;
      /** \brief Number of accepted Delta decays */
      Int_t nDecays;
      /** \brief Number of two-body collisions blocked by Pauli or CDPP */
      Int_t nBlockedCollisions;
      /** \brief Number of decays blocked by Pauli or CDPP */
      Int_t nBlockedDecays;
      /** \brief Number of reflection avatars */
      /** \brief Effective (Coulomb-distorted) impact parameter [fm] */
      Float_t effectiveImpactParameter;

      /// \brief Event involved deltas in the nucleus at the end of the cascade
      Bool_t deltasInside;
      /// \brief Event involved forced delta decays inside the nucleus
      Bool_t forcedDeltasInside;
      /// \brief Event involved forced delta decays outside the nucleus
      Bool_t forcedDeltasOutside;
      /// \brief Event involved cluster decay
      Bool_t clusterDecay;

      /** \brief Time of the first collision [fm/c] */
      Float_t firstCollisionTime;
      /** \brief Cross section of the first collision (mb) */
      Float_t firstCollisionXSec;

      Int_t nReflectionAvatars;
      /** \brief Number of collision avatars */
      Int_t nCollisionAvatars;
      /** \brief Number of decay avatars */
      Int_t nDecayAvatars;

      /// \brief Number of dynamical spectators that were merged back into the projectile remnant
      Int_t nUnmergedSpectators;

      /** \brief Maximum array size for remnants */
      static const Short_t maxSizeRemnants = 10;
      /** \brief Remnant mass number */
      Short_t ARem[maxSizeRemnants];
      /** \brief Remnant charge number */
      Short_t ZRem[maxSizeRemnants];
      /** \brief Remnant excitation energy [MeV] */
      Float_t EStarRem[maxSizeRemnants];
      /** \brief Remnant spin [\f$\hbar\f$] */
      Float_t JRem[maxSizeRemnants];
      /** \brief Remnant kinetic energy [MeV] */
      Float_t EKinRem[maxSizeRemnants];
      /** \brief Remnant momentum, x component [MeV/c] */
      Float_t pxRem[maxSizeRemnants];
      /** \brief Remnant momentum, y component [MeV/c] */
      Float_t pyRem[maxSizeRemnants];
      /** \brief Remnant momentum, z component [MeV/c] */
      Float_t pzRem[maxSizeRemnants];
      /** \brief Remnant momentum polar angle [radians] */
      Float_t thetaRem[maxSizeRemnants];
      /** \brief Remnant momentum azimuthal angle [radians] */
      Float_t phiRem[maxSizeRemnants];
      /** \brief Remnant angular momentum, x component [hbar] */
      Float_t jxRem[maxSizeRemnants];
      /** \brief Remnant angular momentum, y component [hbar] */
      Float_t jyRem[maxSizeRemnants];
      /** \brief Remnant angular momentum, z component [hbar] */
      Float_t jzRem[maxSizeRemnants];

      /** \brief Maximum array size for emitted particles */
      static const Short_t maxSizeParticles = 1000;
      /** \brief Particle mass number */
      Short_t A[maxSizeParticles];
      /** \brief Particle charge number */
      Short_t Z[maxSizeParticles];
      /** \brief Emission time [fm/c] */
      Float_t emissionTime[maxSizeParticles];
      /** \brief Particle kinetic energy [MeV] */
      Float_t EKin[maxSizeParticles];
      /** \brief Particle momentum, x component [MeV/c] */
      Float_t px[maxSizeParticles];
      /** \brief Particle momentum, y component [MeV/c] */
      Float_t py[maxSizeParticles];
      /** \brief Particle momentum, z component [MeV/c] */
      Float_t pz[maxSizeParticles];
      /** \brief Particle momentum polar angle [radians] */
      Float_t theta[maxSizeParticles];
      /** \brief Particle momentum azimuthal angle [radians] */
      Float_t phi[maxSizeParticles];
      /** \brief Origin of the particle
       *
       * Should be -1 for cascade particles, or the number of the remnant for
       * de-excitation particles.
       *
       */
      Short_t origin[maxSizeParticles];
      /** \brief History of the particle
       *
       * Condensed information about the de-excitation chain of a particle. For
       * cascade particles, it is just an empty string. For particles arising
       * from the de-excitation of a cascade remnant, it is a string of
       * characters. Each character represents one or more identical steps in
       * the de-excitation process. The currently defined possible character
       * values and their meanings are the following:
       *
       * e: evaporation product
       * E: evaporation residue
       * m: multifragmentation
       * a: light partner in asymmetric fission or IMF emission
       * A: heavy partner in asymmetric fission or IMF emission
       * f: light partner in fission
       * F: heavy partner in fission
       * s: saddle-to-scission emission
       * n: non-statistical emission (decay)
       */
      std::vector<std::string> history;

#ifdef INCL_INVERSE_KINEMATICS
      /** \brief Particle kinetic energy, in inverse kinematics [MeV] */
      Float_t EKinPrime[maxSizeParticles];
      /** \brief Particle momentum, z component, in inverse kinematics [MeV/c] */
      Float_t pzPrime[maxSizeParticles];
      /** \brief Particle momentum polar angle, in inverse kinematics [radians] */
      Float_t thetaPrime[maxSizeParticles];
#endif // INCL_INVERSE_KINEMATICS

      /** \brief Reset the EventInfo members */
      void reset() {
        Ap = 0;
        Zp = 0;
        At = 0;
        Zt = 0;
        impactParameter = 0.0;
        effectiveImpactParameter = 0.0;
        stoppingTime = 0.0;
        EBalance = 0.0;
        pLongBalance = 0.0;
        pTransBalance = 0.0;
        nCollisions = 0;
        nBlockedCollisions = 0;
        nDecays = 0;
        nBlockedDecays= 0;
        nDecays = 0;
        nCascadeParticles = 0;
        nRemnants = 0;
        nParticles = 0;
        transparent = true;
        forcedCompoundNucleus = false;
	nucleonAbsorption = false;
	pionAbsorption = false;
        forcedDeltasInside = false;
        forcedDeltasOutside = false;
        deltasInside = false;
        clusterDecay = false;
        nUnmergedSpectators = 0;
      }

#ifdef INCL_INVERSE_KINEMATICS
      void fillInverseKinematics(const Double_t gamma);
#endif // INCL_INVERSE_KINEMATICS
    };
}

#endif /* G4INCLEVENTINFO_HH */
