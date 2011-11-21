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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLCascade_hh
#define G4INCLCascade_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLEventAction.hh"
#include "G4INCLPropagationAction.hh"
#include "G4INCLAvatarAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  class INCL {
    public:
      INCL(Config const * const config);
      INCL(IPropagationModel *aPropagationModel);

      ~INCL();

      void setTarget(G4int A, G4int Z);
      G4bool initializeTarget(G4int A, G4int Z);
      const EventInfo &processEvent(Particle *projectile);

      /** \brief Rescale the energies of the outgoing particles.
       *
       * Allow for the remnant recoil energy by rescaling the energy (and
       * momenta) of the outgoing particles.
       */
      void rescaleOutgoingForRecoil();

      /** \brief Run global conservation checks
       *
       * Check that energy and momentum are correctly conserved. If not, issue a
       * warning.
       *
       * Also feeds the balance variables in theEventInfo.
       */
      void globalConservationChecks();

      /** \brief Stopping criterion for the cascade
       *
       * Returns true if the cascade should continue, and false if any of the
       * stopping criteria is satisfied.
       */
      G4bool continueCascade();

      void finaliseGlobalInfo();
      const GlobalInfo &getGlobalInfo() const { return theGlobalInfo; }

      /** \brief Final calculations before returning the global information */
    G4bool processEvent() { return false; }

    std::string configToString() { return theConfig->echo(); }

    private:
      IPropagationModel *propagationModel;
      G4int theA, theZ;
      G4double maxImpactParameter;
      EventAction *eventAction;
      PropagationAction *propagationAction;
      AvatarAction *avatarAction;
      LoggerSlave *theLoggerSlave;
      Config const * const theConfig;

      EventInfo theEventInfo;
      GlobalInfo theGlobalInfo;
  };
}

#endif
