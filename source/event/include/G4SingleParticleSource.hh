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
// G4SingleParticleSource
//
// Class description:
//
// The Single Particle Source is designed to extend the functionality of the
// G4ParticleGun class. It is designed to allow specification of input
// particles in terms of position, direction (or angular) and energy
// distributions. It is used by the General Particle source class
// and it is derived from G4VPrimaryGenerator.
//
// Note on thread safety:
// G4SingleParticleSource instances can be shared among threads.
// GeneratePrimaryVertex is protected via a mutex because underlying
// generators are not assumed to be thread-safe.
// Note that internal status of this class is assumed to be changed by
// master thread (typically via UI commands)
// Only one thread should use the set-methods here.
// If you use the set methods to set defaults in your
// application take care that only one thread is executing them.
// In addition take care of calling these methods before the run is started
// Do not use these setters during the event loop

// Author: Fan Lei, QinetiQ ltd.
// Customer: ESA/ESTEC
// History:
// - 05/02/2004, Fan Lei - Created.
//     Based on the G4GeneralParticleSource class
// - 06/06/2014, Andrea Dotti
//     Added a mutex to protect access to shared resources (data members)
// --------------------------------------------------------------------
#ifndef G4SingleParticleSource_hh
#define G4SingleParticleSource_hh 1

#include "G4VPrimaryGenerator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"

#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"

class G4SingleParticleSource : public G4VPrimaryGenerator
{
  public:

    G4SingleParticleSource();
      // Constructor: initializes variables and instantiates the 
      // messenger and navigator classes

   ~G4SingleParticleSource() override;
      // Destructor: deletes messenger and prints out run information

    void GeneratePrimaryVertex(G4Event *evt) override;
      // Generate the particles initial parameters

    inline G4SPSPosDistribution* GetPosDist() const { return posGenerator; }
      // Return a pointer to the position distribution generator

    inline G4SPSAngDistribution* GetAngDist() const { return angGenerator; }
      // Return a pointer to the angular distribution generator

    inline G4SPSEneDistribution* GetEneDist() const { return eneGenerator; }
      // Return a pointer to the energy distribution generator

    inline G4SPSRandomGenerator* GetBiasRndm() const { return biasRndm; }
      // Return a pointer to the biased random number generator

    void SetVerbosity(G4int);
      // Set the verbosity level

    void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition);
    inline G4ParticleDefinition* GetParticleDefinition() const
           { return definition; }
      // Get/Set the particle definition of the primary track

    inline void SetParticleCharge(G4double aCharge) { charge = aCharge; }
      // Set the charge state of the primary track

    inline void SetParticlePolarization(const G4ThreeVector& aVal)
           { polarization = aVal; }
    inline const G4ThreeVector& GetParticlePolarization() const
           { return polarization; }
      // Set/Get the polarization state of the primary track

    inline void SetParticleTime(G4double aTime) { time = aTime; }
    inline G4double GetParticleTime() const { return time; }
      // Set/Get the Time

    inline void SetNumberOfParticles(G4int i)
           { NumberOfParticlesToBeGenerated = i; }
    inline G4int GetNumberOfParticles() const
           { return NumberOfParticlesToBeGenerated; }
      // Set/get the number of particles to be generated in the primary track

    inline G4ThreeVector GetParticlePosition() const
           { return ParticleProperties.Get().position; }
    inline G4ThreeVector GetParticleMomentumDirection() const
           { return ParticleProperties.Get().momentum_direction; }
    inline G4double GetParticleEnergy() const
           { return ParticleProperties.Get().energy; }
      // Get the position, direction, and energy of the current particle 

  private:

    G4SPSPosDistribution* posGenerator = nullptr;
    G4SPSAngDistribution* angGenerator = nullptr;
    G4SPSEneDistribution* eneGenerator = nullptr;
    G4SPSRandomGenerator* biasRndm = nullptr;

    // Other particle properties
    // These need to be thread-local because a getter for them exits
    //
    struct part_prop_t
    {
      G4ParticleMomentum momentum_direction;
      G4double energy;
      G4ThreeVector position;
      part_prop_t();
    };

    G4Cache<part_prop_t> ParticleProperties;
    G4int NumberOfParticlesToBeGenerated;
    G4ParticleDefinition* definition = nullptr;
    G4double charge;
    G4double time;
    G4ThreeVector polarization;

    G4int verbosityLevel;
      // Verbosity

    G4Mutex mutex;
      // This can be a shared resource.
      // This mutex is uses in GeneratePrimaryVertex
};

#endif
