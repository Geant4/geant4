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

/** \file G4INCLCoulombNone.cc
 * \brief Placeholder class for no Coulomb distortion.
 *
 * Created on: 14 February 2011
 *     Author: Davide Mancusi
 */

#include "G4INCLCoulombNone.hh"

namespace G4INCL {

  void CoulombNone::bringToSurface(Particle * const p, Nucleus const * const n) const {
      ThreeVector momentumUnitVector = p->getMomentum();
      momentumUnitVector /= momentumUnitVector.mag();

      ThreeVector positionTransverse = p->getTransversePosition();
      const G4double impactParameter = positionTransverse.mag();

      G4double radius = n->getSurfaceRadius(p);
      const G4double radius2 = radius*radius;
      G4double distanceZ2 = radius2 - impactParameter * impactParameter;
      if(distanceZ2 < 0.0) distanceZ2 = 0.0;
      const G4double distanceZ = std::sqrt(distanceZ2);

      ThreeVector position = positionTransverse - momentumUnitVector *
        distanceZ;
      p->setPosition(position);
  }

}
