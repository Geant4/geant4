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

/** \file G4INCLICrossSections.hh
 * \brief Abstract interface for the cross-section classes
 *
 * \date 25nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLICROSSSECTIONS_HH
#define G4INCLICROSSSECTIONS_HH

#include "G4INCLParticle.hh"

namespace G4INCL {
  /// \brief Abstract interface for the cross-section classes
  class ICrossSections {
    public:

      ICrossSections() {}
      virtual ~ICrossSections() {}

      /// \brief Elastic particle-particle cross section
      virtual G4double elastic(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Total (elastic+inelastic) particle-particle cross section
      virtual G4double total(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NDelta->NN
      virtual G4double NDeltaToNN(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NN->NDelta
      virtual G4double NNToNDelta(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NN->xpiN
      virtual G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for piN->NDelta
      virtual G4double piNToDelta(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for piN->piNpi
      virtual G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for PiN->EtaN
      virtual G4double piNToEtaN(Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for PiN->OmegaN
      virtual G4double piNToOmegaN(Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for PiN->EtaPrimeN
      virtual G4double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2) = 0;
	  
			   /// \brief Cross section for EtaN->PiN
			   virtual G4double etaNToPiN(Particle const * const p1, Particle const * const p2) = 0;
	  
			   /// \brief Cross section for EtaN->PiPiN
			   virtual G4double etaNToPiPiN(Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for OmegaN->PiN
      virtual G4double omegaNToPiN(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for OmegaN->PiPiN
      virtual G4double omegaNToPiPiN(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for EtaPrimeN->PiN
      virtual G4double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2) = 0;
	  
		   	/// \brief Cross section for NN->NNEta (inclusive)
		   	virtual G4double NNToNNEta(Particle const * const p1, Particle const * const p2) = 0;
	  
		   	/// \brief Cross section for NN->NNEta (exclusive)
		   	virtual G4double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2) = 0;
			
		   	/// \brief Cross section for NN->NNEtaxPi
			   virtual G4double NNToNNEtaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2) = 0;
	  
		   	/// \brief Cross section for N-Delta-Eta production - NNEta Channel
			   virtual G4double NNToNDeltaEta(Particle const * const p1, Particle const * const p2) = 0;
			   
      /// \brief Cross section for NN->NNEta (inclusive)
      virtual G4double NNToNNOmega(Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for NN->NNEta (exclusive)
      virtual G4double NNToNNOmegaExclu(Particle const * const p1, Particle const * const p2) = 0;
			
      /// \brief Cross section for NN->NNEtaxPi
      virtual G4double NNToNNOmegaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2) = 0;
	  
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual G4double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2) = 0;
      
      
      /// \brief elastic scattering for Nucleon-Strange Particles cross sections
      virtual G4double NYelastic(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbelastic(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKelastic(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Nucleon to Stange particles cross sections
      virtual G4double NNToNLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNSK(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNLKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNSKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNLK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNSK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToNNKKb(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Delta to Stange particles cross sections
      virtual G4double NDeltaToNLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NDeltaToNSK(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2) = 0;
      
      virtual G4double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2) = 0;
         
         /// \brief Nucleon-Pion to Stange particles cross sections
         virtual G4double NpiToLK(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToSK(Particle const * const p1, Particle const * const p2) = 0;
		   virtual G4double p_pimToSzKz(Particle const * const p1, Particle const * const p2) = 0;
		   virtual G4double p_pimToSmKp(Particle const * const p1, Particle const * const p2) = 0;
		   virtual G4double p_pizToSzKp(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToLKpi(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToSKpi(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToLK2pi(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToSK2pi(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToNKKb(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Hyperon cross sections
      virtual G4double NLToNS(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NSToNL(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NSToNS(Particle const * const p1, Particle const * const p2) = 0;
         
         /// \brief Nucleon-Kaon inelastic cross sections
         virtual G4double NKToNK(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NKToNKpi(Particle const * const p1, Particle const * const p2) = 0;
         virtual G4double NKToNK2pi(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-antiKaon inelastic cross sections
      virtual G4double NKbToNKb(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToSpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToLpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToS2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToL2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToNKbpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual G4double NKbToNKb2pi(Particle const * const p1, Particle const * const p2) = 0;


      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution
       */
      virtual G4double calculateNNAngularSlope(G4double energyCM, G4int iso) = 0;

  };
}

#endif
