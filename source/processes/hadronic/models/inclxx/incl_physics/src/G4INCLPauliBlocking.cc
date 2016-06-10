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

#include "G4INCLRandom.hh"
#include "G4INCLPauliBlocking.hh"
#include "G4INCLGlobals.hh"

#include "G4INCLPauliStrict.hh"
#include "G4INCLPauliStandard.hh"
#include "G4INCLPauliStrictStandard.hh"
#include "G4INCLPauliGlobal.hh"
#include "G4INCLCDPP.hh"

namespace G4INCL {

  namespace Pauli {

    namespace {
      G4ThreadLocal IPauli * thePauliBlocker = NULL;
      G4ThreadLocal IPauli * theCDPP = NULL;
    }

    IPauli * getBlocker() { return thePauliBlocker; }

    IPauli * getCDPP() { return theCDPP; }

    void setBlocker(IPauli * const pauliBlocker) {
      thePauliBlocker = pauliBlocker;
    }

    void setCDPP(IPauli * const cdpp) {
      theCDPP = cdpp;
    }

    G4bool isBlocked(ParticleList const &modifiedAndCreated, Nucleus const * const nucleus) {
      G4bool isPauliBlocked = false;
      if(thePauliBlocker != 0) {
        isPauliBlocked = thePauliBlocker->isBlocked(modifiedAndCreated, nucleus);
      }

      return isPauliBlocked;
    }

    G4bool isCDPPBlocked(ParticleList const &created, Nucleus const * const nucleus) {
      G4bool isCDPPBlocked = false;
      if(theCDPP != 0) {
        isCDPPBlocked = theCDPP->isBlocked(created, nucleus);
      }

      return isCDPPBlocked;
    }

    void deleteBlockers() {
      delete thePauliBlocker;
      thePauliBlocker=NULL;
      delete theCDPP;
      theCDPP=NULL;
    }

    void initialize(Config const * const aConfig) {
      // Select the Pauli blocking algorithm:
      PauliType pauli = aConfig->getPauliType();
      if(pauli == StrictStatisticalPauli)
        setBlocker(new PauliStrictStandard);
      else if(pauli == StatisticalPauli)
        setBlocker(new PauliStandard);
      else if(pauli == StrictPauli)
        setBlocker(new PauliStrict);
      else if(pauli == GlobalPauli)
        setBlocker(new PauliGlobal);
      else if(pauli == NoPauli)
        setBlocker(NULL);

      if(aConfig->getCDPP())
        setCDPP(new CDPP);
      else
        setCDPP(NULL);

    }
  }
}
