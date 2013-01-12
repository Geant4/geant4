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

/** \file G4INCLInverseInterpolationTable.hh
 * \brief Simple interpolation table for the inverse of a IFunction1D functor
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLINVERSEINTERPOLATIONTABLE_HH_
#define G4INCLINVERSEINTERPOLATIONTABLE_HH_

#include "G4INCLIFunction1D.hh"
#include <algorithm>
#include <functional>
#include <sstream>

namespace G4INCL {

  // Forward declaration
  class InverseInterpolationTable;

  /// \brief Interpolation node
  class InterpolationNode {
    public:
      InterpolationNode(const G4double x0, const G4double y0, const G4double yPrime0) :
        x(x0),
        y(y0),
        yPrime(yPrime0)
    {}

      virtual ~InterpolationNode() {}

      G4bool operator<(const InterpolationNode &rhs) const {
        return (x < rhs.x);
      }

      G4bool operator<=(const InterpolationNode &rhs) const {
        return (x <= rhs.x);
      }

      G4bool operator>(const InterpolationNode &rhs) const {
        return (x > rhs.x);
      }

      G4bool operator>=(const InterpolationNode &rhs) const {
        return (x >= rhs.x);
      }

      /// \brief Overloaded comparison operator for STL algorithms
      friend G4bool operator<(const InterpolationNode &lhs, const G4double rhs) {
        return lhs.x < rhs;
      }

      G4double getX() const { return x; }
      G4double getY() const { return y; }
      G4double getYPrime() const { return yPrime; }

      void setX(const G4double x0) { x=x0; }
      void setY(const G4double y0) { y=y0; }
      void setYPrime(const G4double yPrime0) { yPrime=yPrime0; }

      std::string print() const {
        std::stringstream message;
        message << "x, y, yPrime: " << x << '\t' << y << '\t' << yPrime << std::endl;
        return message.str();
      }

    protected:
      /// \brief abscissa
      G4double x;
      /// \brief function value
      G4double y;
      /// \brief function derivative
      G4double yPrime;
  };

  /// \brief Class for interpolating the inverse of a 1-dimensional function
  class InverseInterpolationTable : public IFunction1D {
    public:
      InverseInterpolationTable(IFunction1D const &f, const unsigned int nNodes=30);
      InverseInterpolationTable(std::vector<G4double> const &x, std::vector<G4double> const &y);
      virtual ~InverseInterpolationTable() {}

      unsigned int getNumberOfNodes() const { return nodes.size(); }

      std::vector<G4double> getNodeAbscissae() const {
        std::vector<G4double> x(nodes.size());
        std::transform(nodes.begin(), nodes.end(), x.begin(),
            std::mem_fun_ref(&InterpolationNode::getX));
        return x;
      }

      std::vector<G4double> getNodeValues() const {
        std::vector<G4double> y(nodes.size());
        std::transform(nodes.begin(), nodes.end(), y.begin(),
            std::mem_fun_ref(&InterpolationNode::getY));
        return y;
      }

      G4double operator()(const G4double x) const;

      std::string print() const;

    private:
      /// \brief Initialise the values of the node derivatives
      void initDerivatives();

      /// \brief Set the function domain
      void setFunctionDomain();

      /// \brief Interpolating nodes
      std::vector<InterpolationNode> nodes;

  };

}

#endif // G4INCLINVERSEINTERPOLATIONTABLE_HH_
