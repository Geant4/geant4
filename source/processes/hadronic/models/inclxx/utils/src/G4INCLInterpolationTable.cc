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

/** \file G4INCLInterpolationTable.cc
 * \brief Simple interpolation table
 *
 * \date 30 January 2014
 * \author Davide Mancusi
 */

// #include <cassert>
#include <algorithm>
#include <functional>
#include "G4INCLInterpolationTable.hh"

namespace G4INCL {

  InterpolationTable::InterpolationTable() : IFunction1D() {}

  InterpolationTable::InterpolationTable(std::vector<G4double> const &x, std::vector<G4double> const &y)
    : IFunction1D(x.front(), x.back())
  {
// assert(x.size()==y.size());
    // Assert that the x vector is sorted
// assert(std::adjacent_find(nodes.begin(), nodes.end(), std::greater<InterpolationNode>()) == nodes.end());

    for(unsigned i = 0; i < x.size(); ++i)
      nodes.push_back(InterpolationNode(x.at(i), y.at(i), 0.));

    initDerivatives();
  }

  std::vector<G4double> InterpolationTable::getNodeAbscissae() const {
    std::vector<G4double> x(nodes.size());
    std::transform(nodes.begin(), nodes.end(), x.begin(),
                   std::mem_fun_ref(&InterpolationNode::getX));
    return x;
  }

  std::vector<G4double> InterpolationTable::getNodeValues() const {
    std::vector<G4double> y(nodes.size());
    std::transform(nodes.begin(), nodes.end(), y.begin(),
                   std::mem_fun_ref(&InterpolationNode::getY));
    return y;
  }

  void InterpolationTable::initDerivatives() {
    for(unsigned i = 0; i < nodes.size()-1; i++) {
      if((nodes.at(i+1).getX() - nodes.at(i).getX()) == 0.0) // Safeguard against division by zero
        nodes[i].setYPrime(0.0);
      else
        nodes[i].setYPrime((nodes.at(i+1).getY() - nodes.at(i).getY())/(nodes.at(i+1).getX() - nodes.at(i).getX()));
    }
    nodes.back().setYPrime(nodes.at(nodes.size()-2).getYPrime()); // Duplicate the last value
  }

  G4double InterpolationTable::operator()(const G4double x) const {
    // Find the relevant interpolation bin
    InterpolationNode xNode(x,0.,0.);
    std::vector<InterpolationNode>::const_iterator iter =
      std::lower_bound(nodes.begin(), nodes.end(), xNode);

    if(iter==nodes.begin())
      return nodes.front().getY();

    if(iter==nodes.end())
      return nodes.back().getY();

    std::vector<InterpolationNode>::const_iterator previousIter = iter - 1;
    const G4double dx = x - previousIter->getX();
    return previousIter->getY() + previousIter->getYPrime()*dx;
  }

  std::string InterpolationTable::print() const {
    std::string message;
    for(std::vector<InterpolationNode>::const_iterator n=nodes.begin(), e=nodes.end(); n!=e; ++n)
      message += n->print();
    return message;
  }

}
