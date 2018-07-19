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

/** \file G4INCLCrossSectionsStrangeness.hh
 * \brief Multipion, mesonic Resonances and strange cross sections
 *
 * \date 1st March 2016
 * \author Jason Hirtz
 */

#ifndef G4INCLCROSSSECTIONSSTRANGENESS_HH
#define G4INCLCROSSSECTIONSSTRANGENESS_HH

#include "G4INCLCrossSectionsMultiPionsAndResonances.hh"
#include "G4INCLConfig.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion, mesonic Resonances and strange cross sections

  class CrossSectionsStrangeness : public CrossSectionsMultiPionsAndResonances {
    public:
      CrossSectionsStrangeness();
	  
      /// \brief second new total particle-particle cross section
      virtual G4double total(Particle const * const p1, Particle const * const p2);
      
      /// \brief second new elastic particle-particle cross section
      virtual G4double elastic(Particle const * const p1, Particle const * const p2);
      
      /// \brief correction to old cross section
      virtual G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2);
      virtual G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2);
	  
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
		G4double p_pimToLK0(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSK(Particle const * const p1, Particle const * const p2);
		G4double p_pipToSpKp(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
		virtual G4double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual G4double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Hyperon quasi-elastic cross sections
      virtual G4double NLToNS(Particle const * const p1, Particle const * const p2);
      virtual G4double NSToNL(Particle const * const p1, Particle const * const p2);
      virtual G4double NSToNS(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Kaon cross sections
      virtual G4double NKToNK(Particle const * const p1, Particle const * const p2);
      virtual G4double NKToNKpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-antiKaon cross sections
      virtual G4double NKbToNKb(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToSpi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToLpi(Particle const * const p1, Particle const * const p2);
		virtual G4double p_kmToL_pz(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      virtual G4double NKbToL2pi(Particle const * const p1, Particle const * const p2);
		virtual G4double p_kmToL_pp_pm(Particle const * const p1, Particle const * const p2);
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
      

  };
}

#endif
