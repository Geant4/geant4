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

/** \file G4INCLCrossSectionsAntiparticles.hh
 * \brief Multipion, mesonic Resonances, strange cross sections and antinucleon as projectile
 *
 * \date 31st March 2023
 * \author Demid Zharenov
 */

#ifndef G4INCLCROSSSECTIONSANTIPARTICLES_HH
#define G4INCLCROSSSECTIONSANTIPARTICLES_HH

#include "G4INCLCrossSectionsStrangeness.hh"
#include "G4INCLConfig.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion, mesonic Resonances and strange cross sections

//  class CrossSectionsAntiparticles : public CrossSectionsMultiPionsAndResonances {
  class CrossSectionsAntiparticles : public CrossSectionsStrangeness {
    public:
      CrossSectionsAntiparticles();
	  
      /// \brief second new total particle-particle cross section
      virtual G4double total(Particle const * const p1, Particle const * const p2);
     
      /// \brief old elastic particle-particle cross section
      virtual G4double elastic(Particle const * const p1, Particle const * const p2);
    
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon cross sections
      virtual G4double NNbarElastic(Particle const* const p1, Particle const* const p2);
      virtual G4double NNbarCEX(Particle const* const p1, Particle const* const p2);

      virtual G4double NNbarToLLbar(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon + pions cross sections
      virtual G4double NNbarToNNbarpi(Particle const* const p1, Particle const* const p2);
      virtual G4double NNbarToNNbar2pi(Particle const* const p1, Particle const* const p2);
      virtual G4double NNbarToNNbar3pi(Particle const* const p1, Particle const* const p2);
     
      /// \brief Nucleon-AntiNucleon total annihilation cross sections
      virtual G4double NNbarToAnnihilation(Particle const* const p1, Particle const* const p2);
         
  protected:
      /// \brief Maximum number of outgoing pions in NN collisions
      static const G4int nMaxPiNN;
	  
      /// \brief Maximum number of outgoing pions in piN collisions
      static const G4int nMaxPiPiN;
	  
      /// \brief Horner coefficients for s11pz
      const HornerC7 s11pzHC;
      /// \brief Horner coefficients for s01pp
      const HornerC8 s01ppHC;
      /// \brief Horner coefficients for s01pz
      const HornerC4 s01pzHC;
      /// \brief Horner coefficients for s11pm
      const HornerC4 s11pmHC;
      /// \brief Horner coefficients for s12pm
      const HornerC5 s12pmHC;
      /// \brief Horner coefficients for s12pp
      const HornerC3 s12ppHC;
      /// \brief Horner coefficients for s12zz
      const HornerC4 s12zzHC;
      /// \brief Horner coefficients for s02pz
      const HornerC4 s02pzHC;
      /// \brief Horner coefficients for s02pm
      const HornerC6 s02pmHC;
      /// \brief Horner coefficients for s12mz
      const HornerC4 s12mzHC;
      
  };
}

#endif
