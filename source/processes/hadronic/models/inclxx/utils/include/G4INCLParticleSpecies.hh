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

/*
 * G4INCLParticleSpecies.hh
 *
 *  \date Nov 25, 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLPARTICLESPECIES_HH_
#define G4INCLPARTICLESPECIES_HH_

#include "G4INCLParticleType.hh"
#include <string>

namespace G4INCL {

  class ParticleSpecies {
    public:
      /// \brief Convert a string to a particle species
      ParticleSpecies() :
        theType(UnknownParticle),
        theA(0),
        theZ(0),
        theS(0)
    {}
      ParticleSpecies(std::string const &pS);
      ParticleSpecies(ParticleType const t);
      ParticleSpecies(const G4int A, const G4int Z);
      ParticleSpecies(const G4int A, const G4int Z, const G4int S);

      ParticleType theType;
      G4int theA, theZ, theS;
	  
      /** \brief Set a PDG Code (MONTE CARLO PARTICLE NUMBERING)
       *
       * \param 
       * \return integer (identifying number for each particle)
       */
      G4int getPDGCode() const;
	  
    private:
      /** \brief Parse a nuclide name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param pS a normalised string (lowercase)
       */
      void parseNuclide(std::string const &pS);

      /** \brief Parse an element name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param pS a normalised string (lowercase)
       * \return true if the parsing succeeded
       */
      G4bool parseElement(std::string const &pS);

      /** \brief Parse a IUPAC element name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param s a normalised string (lowercase)
       * \return true if the parsing succeeded
       */
      G4bool parseIUPACElement(std::string const &s);
	  
  };

}

#endif /* G4INCLPARTICLESPECIES_HH_ */
