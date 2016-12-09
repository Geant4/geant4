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

/** \file G4INCLCrossSectionsTruncatedMultiPions.cc
 * \brief Truncated multipion cross sections
 *
 * \date 10th December 2014
 * \author Davide Mancusi
 */

#include "G4INCLCrossSectionsTruncatedMultiPions.hh"
// #include <cassert>

namespace G4INCL {

  CrossSectionsTruncatedMultiPions::CrossSectionsTruncatedMultiPions(const G4int nPi) :
    nMaxPi(nPi)
  {}

  G4double CrossSectionsTruncatedMultiPions::elastic(Particle const * const p1, Particle const * const p2) {
//    if(!p1->isPion() && !p2->isPion()) {
     if((p1->isNucleon()||p1->isDelta()) && (p2->isNucleon()||p2->isDelta())){
      return NNElastic(p1, p2);
    } 
//  else if (p1->isNucleon() || p2->isNucleon()) {
     else if ((p1->isNucleon() && p2->isPion()) || (p2->isNucleon() && p1->isPion())){
      G4double pielas = piNTot(p1,p2) - piNIne(p1,p2) - CrossSectionsMultiPions::piNToDelta(p1,p2);
      if (pielas < 0.)
        pielas = 0.;
      return pielas;
    } else
      return 0.0;
  }

  G4double CrossSectionsTruncatedMultiPions::piNToDelta(Particle const * const particle1, Particle const * const particle2) {
    G4double sum = CrossSectionsMultiPions::piNToDelta(particle1, particle2);
    if(nMaxPi<=1) {
      for(G4int i=nMaxPi+1; i<=nMaxPiPiN; ++i)
        sum += CrossSectionsMultiPions::piNToxPiN(i, particle1, particle2);
    }
    return sum;
  }

  G4double CrossSectionsTruncatedMultiPions::NNToxPiNN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
    if(xpi<nMaxPi)
      return CrossSectionsMultiPions::NNToxPiNN(xpi, particle1, particle2);
    else if(xpi==nMaxPi) {
      G4double sum = 0.;
      for(G4int i=xpi; i<=nMaxPiNN; ++i)
        sum += CrossSectionsMultiPions::NNToxPiNN(i, particle1, particle2);
      return sum;
    } else
      return 0.;
  }

  G4double CrossSectionsTruncatedMultiPions::piNToxPiN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
    if(xpi<nMaxPi)
      return CrossSectionsMultiPions::piNToxPiN(xpi, particle1, particle2);
    else if(xpi==nMaxPi) {
      G4double sum = 0.;
      for(G4int i=xpi; i<=nMaxPiPiN; ++i)
        sum += CrossSectionsMultiPions::piNToxPiN(i, particle1, particle2);
      return sum;
    } else
      return 0.;
  }

} // namespace G4INCL

