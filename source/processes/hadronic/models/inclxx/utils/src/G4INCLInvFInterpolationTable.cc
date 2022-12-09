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

/** \file G4INCLInvFInterpolationTable.cc
 * \brief Simple interpolation table for the inverse of a IFunction1D functor
 *
 * \date 17 July 2012
 * \author Davide Mancusi
 */

// #include <cassert>
#include <algorithm>
#include <functional>
#include "G4INCLInvFInterpolationTable.hh"

namespace G4INCL {

  InvFInterpolationTable::InvFInterpolationTable(IFunction1D const &f, const unsigned int nNodes) {
// assert(nNodes>2);

    const G4double x0 = f.getXMinimum();
    const G4double x1 = f.getXMaximum();

    // Build the nodes
    G4double last = f(x0);
    InterpolationNode firstNode(last, x0, 0.);
    nodes.push_back(firstNode);
//  G4int skippedNodes = 0;
    for(unsigned i = 1; i < nNodes; i++) {
      const G4double xi = x0 + i*(x1-x0)/((G4double)(nNodes-1));
      // Make sure that the x vector is sorted (corresponding to a monotonous
      // function)
      const G4double value = f(xi);
      if(value <= last) {
//      ++skippedNodes;
        continue;
      }
      InterpolationNode node(value, xi, 0.);
      nodes.push_back(node);
      last = value;
    }

// assert(nNodes==nodes.size()+skippedNodes);

    // Initialise the "derivative" values
    initDerivatives();
  }

}
