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

#include "G4INCLCrossSections.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLLogger.hh"
#include "G4INCLCrossSectionsINCL46.hh"
#include "G4INCLCrossSectionsMultiPions.hh"
#include "G4INCLCrossSectionsTruncatedMultiPions.hh"
#include "G4INCLCrossSectionsMultiPionsAndResonances.hh"
#include "G4INCLCrossSectionsStrangeness.hh"
// #include <cassert>

namespace G4INCL {

  namespace {
    G4ThreadLocal ICrossSections *theCrossSections;
  }

  namespace CrossSections {
    G4double elastic(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->elastic(p1,p2);
    }

    G4double total(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->total(p1,p2);
    }

    G4double NDeltaToNN(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NDeltaToNN(p1,p2);
    }

    G4double NNToNDelta(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNDelta(p1,p2);
    }

	G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2) {
		return theCrossSections->NNToxPiNN(xpi,p1,p2);
	}

	G4double piNToDelta(Particle const * const p1, Particle const * const p2) {
          return theCrossSections->piNToDelta(p1,p2);
	}

	G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2) {
          return theCrossSections->piNToxPiN(xpi,p1,p2);
	}
	  
	G4double piNToEtaN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->piNToEtaN(p1,p2);
	}
	  
	G4double piNToOmegaN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->piNToOmegaN(p1,p2);
	}
	  
	G4double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->piNToEtaPrimeN(p1,p2);
	}
	  
	G4double etaNToPiN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->etaNToPiN(p1,p2);
	}
	  
	G4double etaNToPiPiN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->etaNToPiPiN(p1,p2);
	}
	  
	G4double omegaNToPiN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->omegaNToPiN(p1,p2);
	}
	  
	G4double omegaNToPiPiN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->omegaNToPiPiN(p1,p2);
	}
	  
	G4double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->etaPrimeNToPiN(p1,p2);
	}
	  
	G4double NNToNNEta(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNNEta(p1,p2);
	}
	  
	G4double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNNEtaExclu(p1,p2);
	}
			
	G4double NNToNNEtaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2) {
				return theCrossSections->NNToNNEtaxPi(xpi,p1,p2);
	}
	  
	G4double NNToNDeltaEta(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNDeltaEta(p1,p2);
	}

	  
   G4double NNToNNOmega(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNNOmega(p1,p2);
   }
	  
   G4double NNToNNOmegaExclu(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNNOmegaExclu(p1,p2);
   }
			
   G4double NNToNNOmegaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2) {
				return theCrossSections->NNToNNOmegaxPi(xpi,p1,p2);
   }
	  
   G4double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2) {
		  return theCrossSections->NNToNDeltaOmega(p1,p2);
   }

	
    G4double NYelastic(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NYelastic(p1,p2);
    }
    
    G4double NKbelastic(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbelastic(p1,p2);
    }
    
    G4double NKelastic(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKelastic(p1,p2);
    }
	
    G4double NNToNLK(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNLK(p1,p2);
    }

    G4double NNToNSK(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNSK(p1,p2);
    }

    G4double NNToNLKpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNLKpi(p1,p2);
    }

    G4double NNToNSKpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNSKpi(p1,p2);
    }

    G4double NNToNLK2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNLK2pi(p1,p2);
    }

    G4double NNToNSK2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNSK2pi(p1,p2);
    }

    G4double NNToNNKKb(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNNKKb(p1,p2);
    }
    
    G4double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToMissingStrangeness(p1,p2);
    }

    G4double NDeltaToNLK(Particle const * const p1, Particle const * const p2) {
        return theCrossSections->NDeltaToNLK(p1,p2);
    }
    G4double NDeltaToNSK(Particle const * const p1, Particle const * const p2) {
        return theCrossSections->NDeltaToNSK(p1,p2);
    }
    G4double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2) {
        return theCrossSections->NDeltaToDeltaLK(p1,p2);
    }
    G4double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2) {
        return theCrossSections->NDeltaToDeltaSK(p1,p2);
    }
    
    G4double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2) {
        return theCrossSections->NDeltaToNNKKb(p1,p2);
    }

    G4double NpiToLK(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToLK(p1,p2);
    }

    G4double NpiToSK(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToSK(p1,p2);
    }
    
    G4double p_pimToSmKp(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->p_pimToSmKp(p1,p2);
    }
    
    G4double p_pimToSzKz(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->p_pimToSzKz(p1,p2);
    }
    
	G4double p_pizToSzKp(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->p_pizToSzKp(p1,p2);
    }
	
    G4double NpiToLKpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToLKpi(p1,p2);
    }

    G4double NpiToSKpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToSKpi(p1,p2);
    }

    G4double NpiToLK2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToLK2pi(p1,p2);
    }

    G4double NpiToSK2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToSK2pi(p1,p2);
    }

    G4double NpiToNKKb(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToNKKb(p1,p2);
    }
    
    G4double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NpiToMissingStrangeness(p1,p2);
    }

    G4double NLToNS(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NLToNS(p1,p2);
    }

    G4double NSToNL(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NSToNL(p1,p2);
    }
    
    G4double NSToNS(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NSToNS(p1,p2);
    }

    G4double NKToNK(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKToNK(p1,p2);
    }

    G4double NKToNKpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKToNKpi(p1,p2);
    }

    G4double NKToNK2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKToNK2pi(p1,p2);
    }

    G4double NKbToNKb(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToNKb(p1,p2);
    }

    G4double NKbToSpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToSpi(p1,p2);
    }

    G4double NKbToLpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToLpi(p1,p2);
    }

    G4double NKbToS2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToS2pi(p1,p2);
    }

    G4double NKbToL2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToL2pi(p1,p2);
    }

    G4double NKbToNKbpi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToNKbpi(p1,p2);
    }

    G4double NKbToNKb2pi(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NKbToNKb2pi(p1,p2);
    }


    G4double calculateNNAngularSlope(G4double energyCM, G4int iso) {
      return theCrossSections->calculateNNAngularSlope(energyCM, iso);
    }

    G4double interactionDistancePiN(const G4double projectileKineticEnergy) {
      ThreeVector nullVector;
      ThreeVector unitVector(0., 0., 1.);

      Particle piPlusProjectile(PiPlus, unitVector, nullVector);
      piPlusProjectile.setEnergy(piPlusProjectile.getMass()+projectileKineticEnergy);
      piPlusProjectile.adjustMomentumFromEnergy();
      Particle piZeroProjectile(PiZero, unitVector, nullVector);
      piZeroProjectile.setEnergy(piZeroProjectile.getMass()+projectileKineticEnergy);
      piZeroProjectile.adjustMomentumFromEnergy();
      Particle piMinusProjectile(PiMinus, unitVector, nullVector);
      piMinusProjectile.setEnergy(piMinusProjectile.getMass()+projectileKineticEnergy);
      piMinusProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmapipp = total(&piPlusProjectile, &protonTarget);
      const G4double sigmapipn = total(&piPlusProjectile, &neutronTarget);
      const G4double sigmapi0p = total(&piZeroProjectile, &protonTarget);
      const G4double sigmapi0n = total(&piZeroProjectile, &neutronTarget);
      const G4double sigmapimp = total(&piMinusProjectile, &protonTarget);
      const G4double sigmapimn = total(&piMinusProjectile, &neutronTarget);
      /* We compute the interaction distance from the largest of the pi-N cross
       * sections. Note that this is different from INCL4.6, which just takes the
       * average of the six, and will in general lead to a different geometrical
       * cross section.
       */
      const G4double largestSigma = std::max(sigmapipp, std::max(sigmapipn, std::max(sigmapi0p, std::max(sigmapi0n, std::max(sigmapimp,sigmapimn)))));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    G4double interactionDistanceNN(const ParticleSpecies &aSpecies, const G4double kineticEnergy) {
// assert(aSpecies.theType==Proton || aSpecies.theType==Neutron || aSpecies.theType==Composite);
// assert(aSpecies.theA>0);
      ThreeVector nullVector;
      ThreeVector unitVector(0.,0.,1.);

      const G4double kineticEnergyPerNucleon = kineticEnergy / aSpecies.theA;

      Particle protonProjectile(Proton, unitVector, nullVector);
      protonProjectile.setEnergy(protonProjectile.getMass()+kineticEnergyPerNucleon);
      protonProjectile.adjustMomentumFromEnergy();
      Particle neutronProjectile(Neutron, unitVector, nullVector);
      neutronProjectile.setEnergy(neutronProjectile.getMass()+kineticEnergyPerNucleon);
      neutronProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmapp = total(&protonProjectile, &protonTarget);
      const G4double sigmapn = total(&protonProjectile, &neutronTarget);
      const G4double sigmann = total(&neutronProjectile, &neutronTarget);
      /* We compute the interaction distance from the largest of the NN cross
       * sections. Note that this is different from INCL4.6, which just takes the
       * average of the four, and will in general lead to a different geometrical
       * cross section.
       */
      const G4double largestSigma = std::max(sigmapp, std::max(sigmapn, sigmann));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    G4double interactionDistanceKN(const G4double kineticEnergy) {
      ThreeVector nullVector;
      ThreeVector unitVector(0.,0.,1.);

      Particle kpProjectile(KPlus, unitVector, nullVector);
      kpProjectile.setEnergy(kpProjectile.getMass()+kineticEnergy);
      kpProjectile.adjustMomentumFromEnergy();
      Particle kzProjectile(KZero, unitVector, nullVector);
      kzProjectile.setEnergy(kzProjectile.getMass()+kineticEnergy);
      kzProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmakpp = total(&kpProjectile, &protonTarget);
      const G4double sigmakpn = total(&kpProjectile, &neutronTarget);
      const G4double sigmakzp = total(&kzProjectile, &protonTarget);
      const G4double sigmakzn = total(&kzProjectile, &neutronTarget);
      
      const G4double largestSigma = std::max(sigmakpp, std::max(sigmakpn, std::max(sigmakzp, sigmakzn)));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    G4double interactionDistanceKbarN(const G4double kineticEnergy) {
      ThreeVector nullVector;
      ThreeVector unitVector(0.,0.,1.);

      Particle kmProjectile(KMinus, unitVector, nullVector);
      kmProjectile.setEnergy(kmProjectile.getMass()+kineticEnergy);
      kmProjectile.adjustMomentumFromEnergy();
      Particle kzProjectile(KZeroBar, unitVector, nullVector);
      kzProjectile.setEnergy(kzProjectile.getMass()+kineticEnergy);
      kzProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmakmp = total(&kmProjectile, &protonTarget);
      const G4double sigmakmn = total(&kmProjectile, &neutronTarget);
      const G4double sigmakzp = total(&kzProjectile, &protonTarget);
      const G4double sigmakzn = total(&kzProjectile, &neutronTarget);
      
      const G4double largestSigma = std::max(sigmakmp, std::max(sigmakmn, std::max(sigmakzp, sigmakzn)));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    G4double interactionDistanceYN(const G4double kineticEnergy) {
      ThreeVector nullVector;
      ThreeVector unitVector(0.,0.,1.);

      Particle lProjectile(Lambda, unitVector, nullVector);
      lProjectile.setEnergy(lProjectile.getMass()+kineticEnergy);
      lProjectile.adjustMomentumFromEnergy();
      Particle spProjectile(SigmaPlus, unitVector, nullVector);
      spProjectile.setEnergy(spProjectile.getMass()+kineticEnergy);
      spProjectile.adjustMomentumFromEnergy();
      Particle szProjectile(SigmaZero, unitVector, nullVector);
      szProjectile.setEnergy(szProjectile.getMass()+kineticEnergy);
      szProjectile.adjustMomentumFromEnergy();
      Particle smProjectile(SigmaMinus, unitVector, nullVector);
      smProjectile.setEnergy(smProjectile.getMass()+kineticEnergy);
      smProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmalp = total(&lProjectile, &protonTarget);
      const G4double sigmaln = total(&lProjectile, &neutronTarget);
      const G4double sigmaspp = total(&spProjectile, &protonTarget);
      const G4double sigmaspn = total(&spProjectile, &neutronTarget);
      const G4double sigmaszp = total(&szProjectile, &protonTarget);
      const G4double sigmaszn = total(&szProjectile, &neutronTarget);
      const G4double sigmasmp = total(&smProjectile, &protonTarget);
      const G4double sigmasmn = total(&smProjectile, &neutronTarget);
      
      const G4double largestSigma = std::max(sigmalp, std::max(sigmaln, std::max(sigmaspp, std::max(sigmaspn, std::max(sigmaszp, std::max(sigmaszn, std::max(sigmasmp, sigmasmn)))))));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    void setCrossSections(ICrossSections *c) {
      theCrossSections = c;
    }

    void deleteCrossSections() {
      delete theCrossSections;
      theCrossSections = NULL;
    }

	  void initialize(Config const * const theConfig) {
		  CrossSectionsType crossSections = theConfig->getCrossSectionsType();
		  if(crossSections == INCL46CrossSections)
			  setCrossSections(new CrossSectionsINCL46);
		  else if(crossSections == MultiPionsCrossSections)
			  setCrossSections(new CrossSectionsMultiPions);
		  else if(crossSections == TruncatedMultiPionsCrossSections) {
			  const G4int nMaxPi = theConfig->getMaxNumberMultipions();
			  if(nMaxPi>0)
				  setCrossSections(new CrossSectionsTruncatedMultiPions(nMaxPi));
			  else {
				  INCL_WARN("Truncated multipion cross sections were requested, but the specified maximum\n"
							<< "number of pions is <=0. Falling back to standard multipion cross-sections.\n");
				  setCrossSections(new CrossSectionsMultiPions);
			  }
		  } else if(crossSections == MultiPionsAndResonancesCrossSections)
			  setCrossSections(new CrossSectionsMultiPionsAndResonances);
			else if(crossSections == StrangenessCrossSections)
			  setCrossSections(new CrossSectionsStrangeness);
	  }
  }
}
