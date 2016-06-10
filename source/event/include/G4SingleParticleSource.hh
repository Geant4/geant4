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

///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4SingleParticleSource.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//      Based on the G4GeneralParticleSource class in Geant4 v6.0
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// The Single Particle Source is designed to extend the functionality of the
// G4ParticleGun class. It is designed to allow specification of input
// particles in terms of position, direction (or angular) and energy
// distributions.  It is used by the General Particle source class
// and it is derived from G4VPrimaryGenerator.
//
// Note on thread safety:
//    G4SingleParticleSource instances can be shared among threads. GeneratePrimaryVertex
//    is protected via a mutex because underlying generators are not assumed to be thread-safe.
//    Note that internal status of this class is assumed to be changed by master thread (typically
//    via UI commands).
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4SingleParticleSource ()
//    Constructor: Initializes variables and instantiates the 
//                 Messenger and Navigator classes
//
// ~G4SingleParticleSource ()
//    Destructor:  deletes Messenger and prints out run information.
//
// void GeneratePrimaryVertex(G4Event *evt)
//    Generate the particles initial parameters.
//
//  G4SPSPosDistribution* GetPosDist()
//    Return a pointer to the position distribution generator
//
//  G4SPSAngDistribution* GetAngDist() 
//    Return a pointer to the angular distribution generator
//
//  G4SPSEneDistribution* GetEneDist() 
//     Return a pointer to the energy distribution generator
//
//  G4SPSRandomGenerator* GetBiasRndm() {return biasRndm;};
//     Return a pointer to the biased random number generator
//
//  void SetVerbosity(G4int);
//     Set the verbosity level.
//
//  void SetParticleDefinition ();
//  G4ParticleDefinition * GetParticleDefinition () 
//     Get/Set the particle definition of the primary track
//
//  void SetParticleCharge(G4double aCharge) 
//     set the charge state of the primary track
//
//  inline void SetParticlePolarization (G4ThreeVector aVal) 
//  inline G4ThreeVector GetParticlePolarization ()
//     Set/Get the polarization state of the primary track
//
//  inline void SetParticleTime(G4double aTime)  { particle_time = aTime; };
//  inline G4double GetParticleTime()  { return particle_time; };
//     Set/Get the Time.
//
//  inline void SetNumberOfParticles(G4int i) 
//  inline G4int GetNumberOfParticles() 
//     set/get the number of particles to be generated in the primary track
//
//  inline G4ThreeVector GetParticlePosition()  
//  inline G4ThreeVector GetParticleMomentumDirection()  
//  inline G4double GetParticleEnergy()  
//     get the position, direction, and energy of the current particle 
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef G4SingleParticleSource_h
#define G4SingleParticleSource_h 1

#include "G4VPrimaryGenerator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
//
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"

/** Andrea Dotti Feb 2015
 * Important: This is a shared class between threads.
 * Only one thread should use the set-methods here.
 * Note that this is exactly what is achieved using UI commands.
 * If you use the set methods to set defaults in your
 * application take care that only one thread is executing them.
 * In addition take care of calling these methods before the run is started
 * Do not use these setters during the event loop
 */

class G4SingleParticleSource: public G4VPrimaryGenerator {
public:
	G4SingleParticleSource();
	~G4SingleParticleSource();

	void GeneratePrimaryVertex(G4Event *evt);
	//

	G4SPSPosDistribution* GetPosDist() const {
		return posGenerator;
	}
	;
	G4SPSAngDistribution* GetAngDist() const {
		return angGenerator;
	}
	;
	G4SPSEneDistribution* GetEneDist() const {
		return eneGenerator;
	}
	;
	G4SPSRandomGenerator* GetBiasRndm() const {
		return biasRndm;
	}
	;

	// Set the verbosity level.
	void SetVerbosity(G4int);

	// Set the particle species
	void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
	inline G4ParticleDefinition * GetParticleDefinition() const {
		return definition;
	}
	;

	inline void SetParticleCharge(G4double aCharge) {
	        charge = aCharge;
	}
	;

	// Set polarization
	inline void SetParticlePolarization(G4ThreeVector aVal) {
	  polarization = aVal;
	}
	;
	inline G4ThreeVector GetParticlePolarization() const {
		return polarization;
	}
	;

	// Set Time.
	inline void SetParticleTime(G4double aTime) {
	  time = aTime;
	}
	;
	inline G4double GetParticleTime() const {
		return time;
	}
	;

	inline void SetNumberOfParticles(G4int i) {
	  NumberOfParticlesToBeGenerated = i;
	}
	;
	//
	inline G4int GetNumberOfParticles() const {
		return NumberOfParticlesToBeGenerated;
	}
	;
	inline G4ThreeVector GetParticlePosition() const {
		return ParticleProperties.Get().position;
	}
	;
	inline G4ThreeVector GetParticleMomentumDirection() const {
		return ParticleProperties.Get().momentum_direction;
	}
	;
	inline G4double GetParticleEnergy() const {
		return ParticleProperties.Get().energy;
	}
	;

private:

	G4SPSPosDistribution* posGenerator;
	G4SPSAngDistribution* angGenerator;
	G4SPSEneDistribution* eneGenerator;
	G4SPSRandomGenerator* biasRndm;
	//
	// Other particle properties
	//These need to be thread-local because
	//a getter for them exits
	struct part_prop_t {
	  G4ParticleMomentum momentum_direction; ////////<<<<<<<
	  G4double energy; /////<<<<<
	  G4ThreeVector position; //////////<<<<<<<<<
	  //G4double weight;
	  part_prop_t();
	};
	G4Cache<part_prop_t> ParticleProperties;
        G4int NumberOfParticlesToBeGenerated;
        G4ParticleDefinition * definition;
        G4double charge;
        G4double time;
        G4ThreeVector polarization;

	// Verbosity
	G4int verbosityLevel;

    //This can be a shared resource, this mutex is used in GeneratePrimaryVertex
    G4Mutex mutex;
};

#endif

