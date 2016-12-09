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

#ifndef G4INCLNNOmegaToMultiPionsChannel_hh
#define G4INCLNNOmegaToMultiPionsChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNOmegaToMultiPionsChannel : public IChannel {
    public:
      NNOmegaToMultiPionsChannel(const G4int, Particle *, Particle *);
      virtual ~NNOmegaToMultiPionsChannel();

      void fillFinalState(FinalState *fs);

    private:
      G4int npion;
      G4int iso1; // like isosp, can be changed in isospinRepartition
      G4int iso2; // like isosp, can be changed in isospinRepartition
      G4int isosp[4];
      Particle *particle1, *particle2;

      static const G4double angularSlope;

      void inter2Part(const G4double p);
      void pn_ppPim();
      void pn_pnPi0();
      void pn_nnPip();
      void pp_nnPipPip();
      void nn_ppPimPim();
      void pn_pnPipPim();
      void pn_pnPi0Pi0();
      void pp_ppPipPim();
      void nn_nnPipPim();
      void pp_ppPi0Pi0();
      void nn_nnPi0Pi0();
      void pp_pnPipPi0();
      void pn_ppPimPi0();
      void pn_nnPipPi0();
      void nn_pnPimPi0();
      void pp_pnPipPi0Pi0();
      void nn_pnPimPi0Pi0();
      void pn_nnPipPi0Pi0();
      void pp_ppPipPimPi0();
      void nn_nnPipPimPi0();
      void pp_ppPi0Pi0Pi0();
      void nn_nnPi0Pi0Pi0();
      void pp_pnPipPipPim();
      void pp_nnPipPipPi0();
      void pn_ppPimPi0Pi0();
      void pn_ppPimPimPip();
      void pn_pnPi0PipPim();
      void pn_pnPi0Pi0Pi0();
      void pn_nnPipPipPim();
      void nn_pnPipPimPim();
      void nn_ppPimPimPi0();
      void pp_nnPipPipPi0Pi0();
      void pp_nnPipPipPipPim();
      void nn_ppPi0Pi0PimPim();
      void nn_ppPipPimPimPim();
      void pp_ppPi0Pi0Pi0Pi0();
      void nn_nnPi0Pi0Pi0Pi0();
      void pn_pnPi0Pi0Pi0Pi0();
      void pp_ppPipPi0Pi0Pim();
      void nn_nnPipPi0Pi0Pim();
      void pn_pnPipPi0Pi0Pim();
      void pp_ppPipPipPimPim();
      void nn_nnPipPipPimPim();
      void pn_pnPipPipPimPim();
      void pp_pnPipPi0Pi0Pi0();
      void pn_nnPipPi0Pi0Pi0();
      void pp_nnPipPi0Pi0Pi0();
      void pp_pnPipPipPi0Pim();
      void pn_nnPipPipPi0Pim();
      void pp_nnPipPipPi0Pim();
      void nn_pnPi0Pi0Pi0Pim();
      void pn_ppPi0Pi0Pi0Pim();
      void nn_pnPipPi0PimPim();
      void pn_ppPipPi0PimPim();
      void isospinRepartition();

      INCL_DECLARE_ALLOCATION_POOL(NNOmegaToMultiPionsChannel);
  };
}

#endif
