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

/** \file G4INCLCrossSectionsMultiPions.hh
 * \brief Cross sections used in INCL Multipions
 *
 * \date 26th November 2013
 * \author Jean-Christophe David
 */

#ifndef G4INCLCROSSSECTIONSMULTIPIONS_HH
#define G4INCLCROSSSECTIONSMULTIPIONS_HH

#include "G4INCLICrossSections.hh"
#include "G4INCLHornerFormEvaluator.hh"

namespace G4INCL {
  /// \brief Cross sections used in INCL Multipions

  class CrossSectionsMultiPions : public ICrossSections{
    public:
      CrossSectionsMultiPions();

      /// \brief Elastic particle-particle cross section
      virtual G4double elastic(Particle const * const p1, Particle const * const p2);

      /// \brief Total (elastic+inelastic) particle-particle cross section
      virtual G4double total(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for NDelta->NN
      virtual G4double NDeltaToNN(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - NN Channel
      virtual G4double NNToNDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - piN Channel
      virtual G4double piNToDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - piN Channel
      virtual G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - NN Channel
      virtual G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2);

      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution
       */
      virtual G4double calculateNNAngularSlope(G4double energyCM, G4int iso);
	  
      /// \brief Cross sections for mesonic resonance production - piN Channel
      virtual G4double piNToEtaN(Particle const * const p1, Particle const * const p2);
      virtual G4double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      virtual G4double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
			
			   /// \brief Cross sections for mesonic resonance absorption on nucleon - pipiN Channel
      virtual G4double etaNToPiPiN(Particle const * const p1, Particle const * const p2);			
      virtual G4double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);	  
			
      /// \brief Cross sections for mesonic resonance absorption on nucleon - piN Channel
      virtual G4double etaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual G4double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual G4double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);
	  
	     /// \brief Cross section for Eta production - NN entrance channel
      virtual G4double NNToNNEta(Particle const * const particle1, Particle const * const particle2);
			
		   	/// \brief Cross section for Eta production  (exclusive) - NN entrance channel
		   	virtual G4double NNToNNEtaExclu(Particle const * const particle1, Particle const * const particle2);
	  
		   	/// \brief Cross section for X pion production - NNEta Channel
			   virtual G4double NNToNNEtaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
		   	/// \brief Cross section for N-Delta-Eta production - NNEta Channel
		   	virtual G4double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);			
	  
      /// \brief Cross section for Eta production - NN entrance channel
      virtual G4double NNToNNOmega(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
      virtual G4double NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2);
	  
      /// \brief Cross section for X pion production - NNEta Channel
      virtual G4double NNToNNOmegaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual G4double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
      
      
      /// \brief elastic scattering for Nucleon-Strange Particles cross sections
      virtual G4double NYelastic(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbelastic(Particle const * const p1, Particle const * const p2);
      virtual G4double NKelastic(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Nucleon to Stange particles cross sections
      virtual G4double NNToNLK(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNSK(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNLKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNSKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNLK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNSK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NNToNNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual G4double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Delta to Stange particles cross sections
      virtual G4double NDeltaToNLK(Particle const * const p1, Particle const * const p2);
      virtual G4double NDeltaToNSK(Particle const * const p1, Particle const * const p2);
      virtual G4double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2);
      virtual G4double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2);
      
      virtual G4double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Pion to Stange particles cross sections
      virtual G4double NpiToLK(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSK(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual G4double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Hyperon cross sections
      virtual G4double NLToNS(Particle const * const p1, Particle const * const p2);
      virtual G4double NSToNL(Particle const * const p1, Particle const * const p2);
      virtual G4double NSToNS(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Kaon quasi-elastic and inelastic cross sections
      virtual G4double NKToNK(Particle const * const p1, Particle const * const p2);
      virtual G4double NKToNKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-antiKaon quasi-elastic and inelastic cross sections
      virtual G4double NKbToNKb(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToSpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToLpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToL2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);


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

      /// \brief Internal implementation of the NN elastic cross section
      G4double NNElastic(Particle const * const part1, Particle const * const part2);

      /// \brief Internal implementation of the NN elastic cross section with fixed isospin
      G4double NNElasticFixed(const G4double s, const G4int i);

      /// \brief Internal implementation of the NN total cross section
      G4double NNTot(Particle const * const part1, Particle const * const part2);

      /// \brief Internal implementation of the NN total cross section with fixed isospin
      G4double NNTotFixed(const G4double s, const G4int i);

      /// \brief Internal implementation of the isospin dependent NN reaction cross section
      G4double NNInelasticIso(const G4double ener, const G4int iso);

      /// \brief Cross section for direct 1-pion production + delta production - NN entrance channel
      virtual G4double NNOnePiOrDelta(const G4double ener, const G4int iso, const G4double xsiso);
      /// \brief Cross section for direct 2-pion production - NN entrance channel
      virtual G4double NNTwoPi(const G4double ener, const G4int iso, const G4double xsiso);
      /// \brief Cross section for direct 3-pion production - NN entrance channel
      virtual G4double NNThreePi(const G4double ener, const G4int iso, const G4double xsiso, const G4double xs1pi, const G4double xs2pi);

      /// \brief Cross section for direct 1-pion production - NN entrance channel
      virtual G4double NNOnePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 1-pion production - NN entrance channel
      virtual G4double NNOnePiOrDelta(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 2-pion production - NN entrance channel
      virtual G4double NNTwoPi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 3-pion production - NN entrance channel
      virtual G4double NNThreePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 4-pion production - NN entrance channel
      virtual G4double NNFourPi(Particle const * const part1, Particle const * const part2);

      /// \brief Internal function for pion cross sections
      G4double spnPiPlusPHE(const G4double x);
      /// \brief Internal function for pion cross sections
      G4double spnPiMinusPHE(const G4double x);
      G4double piNIne(Particle const * const p1, Particle const * const p2);
      G4double piNTot(Particle const * const p1, Particle const * const p2);
      G4double piNTopiN(Particle const * const p1, Particle const * const p2);
      G4double piPluspIne(Particle const * const p1, Particle const * const p2);
	  G4double piMinuspIne(Particle const * const p1, Particle const * const p2);
      G4double piPluspOnePi(Particle const * const p1, Particle const * const p2);
	  G4double piMinuspOnePi(Particle const * const p1, Particle const * const p2);
      G4double piPluspTwoPi(Particle const * const p1, Particle const * const p2);
	  G4double piMinuspTwoPi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for One (more) pion production - piN entrance channel
      virtual G4double piNOnePi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Two (more) pion production - piN entrance channel
      virtual G4double piNTwoPi(Particle const * const p1, Particle const * const p2);
	  
	  
  };
}

#endif
