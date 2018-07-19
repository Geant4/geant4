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

/** \file G4INCLCrossSectionsMultiPionsAndResonances.hh
 * \brief Multipion and mesonic Resonances cross sections
 *
 * \date 4th February 2014
 * \author Jean-Christophe David
 */

#ifndef G4INCLCROSSSECTIONSMULTIPIONSANDRESONANCES_HH
#define G4INCLCROSSSECTIONSMULTIPIONSANDRESONANCES_HH

#include "G4INCLCrossSectionsMultiPions.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion and mesonic Resonances cross sections

  class CrossSectionsMultiPionsAndResonances : public CrossSectionsMultiPions {
    public:
      CrossSectionsMultiPionsAndResonances();
	  
      /// \brief new elastic particle-particle cross section
      virtual G4double elastic(Particle const * const p1, Particle const * const p2);
	  
      /// \brief new total particle-particle cross section
      virtual G4double total(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for X pion production - piN Channel (modified due to the mesonic resonances)
      virtual G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance production - piN Channel
      virtual G4double piNToEtaN(Particle const * const p1, Particle const * const p2);
      virtual G4double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      virtual G4double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance absorption on nucleon - piN Channel
      virtual G4double etaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual G4double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual G4double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);

			   /// \brief Cross sections for mesonic resonance absorption on nucleon - pipiN Channel
			   virtual G4double etaNToPiPiN(Particle const * const p1, Particle const * const p2);
			
      /// \brief Cross section for Eta production (inclusive) - NN entrance channel
      virtual G4double NNToNNEta(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
      virtual G4double NNToNNEtaExclu(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Omega production (inclusive) - NN entrance channel
      virtual G4double NNToNNOmega(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Omega production  (exclusive) - NN entrance channel
      virtual G4double NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2);
	  
      /// \brief Cross section for X pion production - NN Channel
      virtual G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
		   	/// \brief Cross section for X pion production - NNEta Channel
		   	virtual G4double NNToNNEtaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
			   /// \brief Cross section for N-Delta-Eta production - NNEta Channel
			   virtual G4double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for X pion production - NNOmega Channel
      virtual G4double NNToNNOmegaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for N-Delta-Eta production - NNOmega Channel
      virtual G4double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
			
	  
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
	  
      /// \brief One over threshold for s11pz
      static const G4double s11pzOOT;
      /// \brief One over threshold for s01pp
      static const G4double s01ppOOT;
      /// \brief One over threshold for s01pz
      static const G4double s01pzOOT;
      /// \brief One over threshold for s11pm
      static const G4double s11pmOOT;
      /// \brief One over threshold for s12pm
      static const G4double s12pmOOT;
      /// \brief One over threshold for s12pp
      static const G4double s12ppOOT;
      /// \brief One over threshold for s12zz
      static const G4double s12zzOOT;
      /// \brief One over threshold for s02pz
      static const G4double s02pzOOT;
      /// \brief One over threshold for s02pm
      static const G4double s02pmOOT;
      /// \brief One over threshold for s12mz
      static const G4double s12mzOOT;
	  
	  
      /// \brief Internal function for pion cross sections
	  G4double piMinuspToEtaN(Particle const * const p1, Particle const * const p2);
	  G4double piMinuspToEtaN(const G4double ECM);
	  G4double piMinuspToOmegaN(Particle const * const p1, Particle const * const p2);
	  G4double piMinuspToOmegaN(const G4double ECM);
//      G4double piPluspOnePi(Particle const * const p1, Particle const * const p2);
//	  G4double piMinuspOnePi(Particle const * const p1, Particle const * const p2);
//      G4double piPluspTwoPi(Particle const * const p1, Particle const * const p2);
//	  G4double piMinuspTwoPi(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for One (more) pion production - piN entrance channel
//      virtual G4double piNOnePi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Two (more) pion production - piN entrance channel
//      virtual G4double piNTwoPi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Three (more) pion production - piN entrance channel
      ///virtual G4double piNThreePi(Particle const * const p1, Particle const * const p2);
	  
		   	/// \brief Isotopic Cross section for Eta production (inclusive) - NN entrance channel
		   	virtual G4double NNToNNEtaIso(const G4double ener, const G4int iso);
	  
		   /// \brief Isotopic Cross section for Eta production (exclusive) - NN entrance channel
		   	virtual G4double NNToNNEtaExcluIso(const G4double ener, const G4int iso);
			
		   	/// \brief Cross section for direct 1-pion production - NNEta  channel
			   virtual G4double NNToNNEtaOnePi(Particle const * const part1, Particle const * const part2);
						/// \brief Cross section for direct 1-pion production - NNEta  channel
		   	virtual G4double NNToNNEtaOnePiOrDelta(Particle const * const part1, Particle const * const part2);
			   /// \brief Cross section for direct 2-pion production - NNEta  channel
		   	virtual G4double NNToNNEtaTwoPi(Particle const * const part1, Particle const * const part2);
   			/// \brief Cross section for direct 3-pion production - NNEta  channel
	   		virtual G4double NNToNNEtaThreePi(Particle const * const part1, Particle const * const part2);
	   		/// \brief Cross section for direct 4-pion production - NNEta  channel
	   		virtual G4double NNToNNEtaFourPi(Particle const * const part1, Particle const * const part2);

	  
      /// \brief Isotopic Cross section for Omega production (inclusive) - NN entrance channel
      virtual G4double NNToNNOmegaIso(const G4double ener, const G4int iso);
	  
      /// \brief Isotopic Cross section for Omega production (exclusive) - NN entrance channel
      virtual G4double NNToNNOmegaExcluIso(const G4double ener, const G4int iso);
			
      /// \brief Cross section for direct 1-pion production - NNOmega  channel
      virtual G4double NNToNNOmegaOnePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 1-pion production - NNOmega  channel
      virtual G4double NNToNNOmegaOnePiOrDelta(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 2-pion production - NNOmega  channel
      virtual G4double NNToNNOmegaTwoPi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 3-pion production - NNOmega  channel
      virtual G4double NNToNNOmegaThreePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 4-pion production - NNOmega  channel
      virtual G4double NNToNNOmegaFourPi(Particle const * const part1, Particle const * const part2);
   
			
		   	/// \brief Cross sections for mesonic resonance absorption on nucleon - elastic Channel
		   	virtual G4double etaNElastic(Particle const * const p1, Particle const * const p2);					
      virtual G4double omegaNElastic(Particle const * const p1, Particle const * const p2);					

			
      /// \brief Cross sections for mesonic resonance absorption on nucleon - inelastic Channel
      virtual G4double omegaNInelastic(Particle const * const p1, Particle const * const p2);					
			
      /// \brief Cross sections for omega-induced 2Pi emission on nucleon
      virtual G4double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);					
   
		};
}

#endif
