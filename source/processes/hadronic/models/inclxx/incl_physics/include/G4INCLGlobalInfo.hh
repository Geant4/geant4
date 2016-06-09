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

/** \file G4INCLGlobalInfo.hh
 * \brief Simple container for output of calculation-wide results.
 *
 * Contains the final results of an INCL calculation.
 *
 * Created on: 21 January 2011
 *     Author: Davide Mancusi
 */

#ifndef G4INCLGLOBALINFO_HH
#define G4INCLGLOBALINFO_HH 1

#ifdef INCL_ROOT_USE
#include <Rtypes.h>
#endif

namespace G4INCL {
#ifndef INCL_ROOT_USE
    typedef G4int Int_t;
    typedef short Short_t;
    typedef G4float Float_t;
#endif

    struct GlobalInfo {
      GlobalInfo() :
        nShots(0), nTransparents(0),
        geometricCrossSection(0.0), reactionCrossSection(0.0),
        Ap(0), Zp(0), At(0), Zt(0), Ep(0.0)
      {};

      /** \brief Number of shots */
      Int_t nShots;
      /** \brief Number of transparent shots */
      Int_t nTransparents;
      /** \brief Geometric cross section */
      Float_t geometricCrossSection;
      /** \brief Calculated reaction cross section */
      Float_t reactionCrossSection;
      /** \brief Error on the calculated reaction cross section */
      Float_t errorReactionCrossSection;

      // TODO: echo all the input parameters here
      /** \brief Projectile mass number given as input */
      Short_t Ap;
      /** \brief Projectile charge number given as input */
      Short_t Zp;
      /** \brief Target mass number given as input */
      Short_t At;
      /** \brief Target charge number given as input */
      Short_t Zt;
      /** \brief Projectile kinetic energy given as input */
      Float_t Ep;

      /** \brief Maximum model name size */
      static const Short_t maxModelSize = 200;
      /** \brief Name of the cascade model */
      char cascadeModel[maxModelSize];
      /** \brief Name of the de-excitation model */
      char deexcitationModel[maxModelSize];
    };
}

#endif /* G4INCLGLOBALINFO_HH */
