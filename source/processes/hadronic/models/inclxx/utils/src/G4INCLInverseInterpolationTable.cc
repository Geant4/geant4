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

/** \file G4INCLInverseInterpolationTable.cc
 * \brief Simple interpolation table for the inverse of a IFunction1D functor
 *
 * \date 17 July 2012
 * \author Davide Mancusi
 */

// #include <cassert>
#include <algorithm>
#include <functional>
#include "G4INCLInverseInterpolationTable.hh"

namespace G4INCL {

  InverseInterpolationTable::InverseInterpolationTable(IFunction1D const &f, const unsigned int nNodes) {
// assert(nNodes>2);

    const G4double x0 = f.getXMinimum();
    const G4double x1 = f.getXMaximum();

    // Build the nodes
    G4double last = f(x0);
    InterpolationNode firstNode(last, x0, 0.);
    nodes.push_back(firstNode);
    G4int skippedNodes = 0;
    for(unsigned i = 1; i < nNodes; i++) {
      const G4double xi = x0 + i*(x1-x0)/((G4double)(nNodes-1));
      // Make sure that the x vector is sorted (corresponding to a monotonous
      // function)
      const G4double value = f(xi);
      if(value <= last) {
        ++skippedNodes;
        continue;
      }
      InterpolationNode node(value, xi, 0.);
      nodes.push_back(node);
      last = value;
    }

// assert(nNodes==nodes.size()+skippedNodes);

    // Initialise the "derivative" values
    initDerivatives();
    setFunctionDomain();
  }

  InverseInterpolationTable::InverseInterpolationTable(std::vector<G4double> const &x, std::vector<G4double> const &y) {
// assert(x.size()==y.size());
    // Assert that the x vector is sorted (corresponding to a monotonous
    // function
// assert(std::adjacent_find(nodes.begin(), nodes.end(), std::greater<InterpolationNode>()) == nodes.end());

    for(unsigned i = 0; i < x.size(); ++i)
      nodes.push_back(InterpolationNode(x.at(i), y.at(i), 0.));

    initDerivatives();
    setFunctionDomain();
  }

  void InverseInterpolationTable::initDerivatives() {
    for(unsigned i = 0; i < nodes.size()-1; i++) {
      if((nodes.at(i+1).getX() - nodes.at(i).getX()) == 0.0) // Safeguard against division by zero
        nodes[i].setYPrime(0.0);
      else
        nodes[i].setYPrime((nodes.at(i+1).getY() - nodes.at(i).getY())/(nodes.at(i+1).getX() - nodes.at(i).getX()));
    }
    nodes.back().setYPrime(nodes.at(nodes.size()-2).getYPrime()); // Duplicate the last value
  }

  void InverseInterpolationTable::setFunctionDomain() {
    // Set the function domain
    if(nodes.front()>nodes.back()) {
      xMin = nodes.back().getX();
      xMax = nodes.front().getX();
    } else {
      xMin = nodes.front().getX();
      xMax = nodes.back().getX();
    }
  }

  G4double InverseInterpolationTable::operator()(const G4double x) const {
    // Find the relevant interpolation bin
    std::vector<InterpolationNode>::const_iterator iter =
      std::lower_bound(nodes.begin(), nodes.end(), x);

    if(iter==nodes.begin())
      return nodes.front().getY();

    if(iter==nodes.end())
      return nodes.back().getY();

    std::vector<InterpolationNode>::const_iterator previousIter = iter - 1;
    const G4double dx = x - previousIter->getX();
    return previousIter->getY() + previousIter->getYPrime()*dx;
  }

  std::string InverseInterpolationTable::print() const {
    std::string message;
    for(std::vector<InterpolationNode>::const_iterator n=nodes.begin(); n!=nodes.end(); ++n)
      message += n->print();
    return message;
  }

}
