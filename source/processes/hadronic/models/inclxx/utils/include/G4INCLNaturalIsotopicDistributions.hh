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

/** \file G4INCLNaturalIsotopicDistributions.hh
 * \brief Classes that stores isotopic abundances
 *
 * \date 21st October 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLISOTOPICDISTRIBUTION_HH_
#define G4INCLISOTOPICDISTRIBUTION_HH_

#include <vector>
#include <map>

namespace G4INCL {

  /// \brief Holds an isotope and an abundance
  struct Isotope {
    Isotope(const G4int A, const G4double abundance);
    G4int theA;
    G4double theAbundance;
  };

  typedef std::vector<Isotope> IsotopeVector;
  typedef IsotopeVector::iterator IsotopeIter;

  /// \brief Class that stores isotopic abundances for a given element
  class IsotopicDistribution {
    public:

      /// \brief Constructor
      IsotopicDistribution(IsotopeVector const &aVector);

      /// \brief Draw a random isotope based on the abundance vector
      G4int drawRandomIsotope() const;

      /// \brief Get the isotope vector
      IsotopeVector const &getIsotopes() const;

    private:
      IsotopeVector theIsotopes;
  };

  /// \brief Class that stores isotopic abundances for a given element
  class NaturalIsotopicDistributions {
    public:

      /// \brief Constructor
      NaturalIsotopicDistributions();

      /** \brief Draw a random isotope
       *
       * \param Z the element number
       */
      G4int drawRandomIsotope(G4int const Z) const;

      /** \brief Get an isotopic distribution
       *
       * \param Z the element number
       */
      IsotopicDistribution const &getIsotopicDistribution(G4int const Z) const;

    private:
      std::map<G4int, IsotopicDistribution> theDistributions;
  };

}
#endif // G4INCLISOTOPICDISTRIBUTION_HH_
