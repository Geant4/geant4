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
// G4DynamicParticle
//
// Class description:
//
// A G4DynamicParticle aggregates the information to describe the dynamics
// of a G4Particle, such as energy, momentum, polarization and proper time,
// as well as the "particle definition", holding all the static information.
// It contains the purely dynamic aspects of a moving particle.

// History:
// - 2 December 1995, G.Cosmo - first design, based on object model.
// - 29 January 1996, M.Asai - first implementation.
// - 1996 - 2007,     H.Kurashige - revisions.
// - 15 March 2019,   M.Novak - log-kinetic energy value is computed only
//                    on demand if its stored value is not up-to-date.
// --------------------------------------------------------------------
#ifndef G4DynamicParticle_hh
#define G4DynamicParticle_hh 1

#include "G4Allocator.hh"
#include "G4ElectronOccupancy.hh"
#include "G4Log.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"  // NOTE: means "momentum direction" not "momentum vector". It is a G4ThreeVector
#include "G4ios.hh"
#include "globals.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include <cmath>

class G4PrimaryParticle;
class G4DecayProducts;

class G4DynamicParticle
{
  public:
    //- constructors
    G4DynamicParticle();

    G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                      const G4ThreeVector& aMomentumDirection, G4double aKineticEnergy);
    G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                      const G4ThreeVector& aParticleMomentum);
    G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                      const G4LorentzVector& aParticleMomentum);
    G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition, G4double aTotalEnergy,
                      const G4ThreeVector& aParticleMomentum);
    G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                      const G4ThreeVector& aMomentumDirection, G4double aKineticEnergy,
                      const G4double dynamicalMass);

    G4DynamicParticle(const G4DynamicParticle& right);

    //- destructor
    ~G4DynamicParticle();

    //- operators
    G4DynamicParticle& operator=(const G4DynamicParticle& right);
    G4bool operator==(const G4DynamicParticle& right) const;
    G4bool operator!=(const G4DynamicParticle& right) const;

    //- Move constructor & operator
    G4DynamicParticle(G4DynamicParticle&& from);
    G4DynamicParticle& operator=(G4DynamicParticle&& from);

    //- new/delete operators are oberloded to use G4Allocator
    inline void* operator new(size_t);
    inline void operator delete(void* aDynamicParticle);

    //- Set/Get methods

    // Returns the normalized direction of the momentum
    inline const G4ThreeVector& GetMomentumDirection() const;

    // Sets the normalized direction of the momentum
    inline void SetMomentumDirection(const G4ThreeVector& aDirection);

    // Sets the normalized direction of the momentum by coordinates
    inline void SetMomentumDirection(G4double px, G4double py, G4double pz);

    // Returns the current particle momentum vector
    inline G4ThreeVector GetMomentum() const;

    // set the current particle momentum vector
    void SetMomentum(const G4ThreeVector& momentum);

    // Returns the current particle energy-momentum 4vector
    inline G4LorentzVector Get4Momentum() const;

    // Set the current particle energy-momentum 4vector
    void Set4Momentum(const G4LorentzVector& momentum);

    // Returns the module of the momentum vector
    inline G4double GetTotalMomentum() const;

    // Returns the total energy of the particle
    inline G4double GetTotalEnergy() const;

    // Returns the kinetic energy of a particle
    inline G4double GetKineticEnergy() const;

    // Returns:
    // - natural logarithm of the particle kinetic energy (E_k) if E_k > 0
    // - LOG_EKIN_MIN otherwise
    inline G4double GetLogKineticEnergy() const;

    // Sets the kinetic energy of a particle
    inline void SetKineticEnergy(G4double aEnergy);

    // Access Lorentz beta
    inline G4double GetBeta() const;

    // Returns the current particle proper time
    inline G4double GetProperTime() const;

    // Set the current particle Proper Time
    inline void SetProperTime(G4double);

    // Set/Get polarization vector
    inline const G4ThreeVector& GetPolarization() const;
    inline void SetPolarization(const G4ThreeVector&);
    inline void SetPolarization(G4double polX, G4double polY, G4double polZ);

    // Set/Get dynamical mass
    // The dynamical mass is set to PDG mass in default
    inline G4double GetMass() const;
    inline void SetMass(G4double mass);

    // Set/Get dynamical charge
    // The dynamical mass is set to PDG charge in default
    inline G4double GetCharge() const;
    inline void SetCharge(G4double charge);
    inline void SetCharge(G4int chargeInUnitOfEplus);

    // Set/Get dynamical spin
    // The dynamical spin is set to PDG spin in default
    inline G4double GetSpin() const;
    inline void SetSpin(G4double spin);
    inline void SetSpin(G4int spinInUnitOfHalfInteger);

    // Set/Get dynamical MagneticMoment
    // The dynamical mass is set to PDG MagneticMoment in default
    inline G4double GetMagneticMoment() const;
    inline void SetMagneticMoment(G4double magneticMoment);

    // Get electron occupancy
    // ElectronOccupancy is valid only if the particle is ion
    inline const G4ElectronOccupancy* GetElectronOccupancy() const;
    inline G4int GetTotalOccupancy() const;
    inline G4int GetOccupancy(G4int orbit) const;
    inline void AddElectron(G4int orbit, G4int number = 1);
    inline void RemoveElectron(G4int orbit, G4int number = 1);

    // Set/Get particle definition
    inline const G4ParticleDefinition* GetParticleDefinition() const;
    void SetDefinition(const G4ParticleDefinition* aParticleDefinition);

    // Following method of GetDefinition() remains
    // because of backward compatiblity. May be removed in future
    inline G4ParticleDefinition* GetDefinition() const;

    // Set/Get pre-assigned decay channel
    inline const G4DecayProducts* GetPreAssignedDecayProducts() const;
    inline void SetPreAssignedDecayProducts(G4DecayProducts* aDecayProducts);

    // Set/Get pre-assigned proper time when the particle will decay
    inline G4double GetPreAssignedDecayProperTime() const;
    inline void SetPreAssignedDecayProperTime(G4double);

    // Print out information
    // - mode 0 : default )(minimum)
    // - mode 1 : 0 + electron occupancy
    void DumpInfo(G4int mode = 0) const;

    // Set/Get controle flag for output message
    // - 0: Silent
    // - 1: Warning message
    // - 2: More
    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;

    inline void SetPrimaryParticle(G4PrimaryParticle* p);
    inline void SetPDGcode(G4int c);

    // Return the pointer to the corresponding G4PrimaryParticle object
    // if this particle is a primary particle OR is defined as a
    // pre-assigned decay product. Otherwise return nullptr.
    inline G4PrimaryParticle* GetPrimaryParticle() const;

    // Return the PDG code of this particle. If the particle is known to
    // Geant4, its PDG code defined in G4ParticleDefinition is returned.
    // If it is unknown (i.e. PDG code in G4ParticleDefinition is 0), the
    // PDG code defined in the corresponding primary particle or
    // pre-assigned decay product will be returned if available.
    // Otherwise (e.g. for geantino) returns 0.
    inline G4int GetPDGcode() const;

  protected:
    void AllocateElectronOccupancy();
    G4double GetElectronMass() const;

  private:
    inline void ComputeBeta() const;

    // The normalized momentum vector
    G4ThreeVector theMomentumDirection;

    G4ThreeVector thePolarization;

    // Contains the static information of this particle
    const G4ParticleDefinition* theParticleDefinition = nullptr;

    G4ElectronOccupancy* theElectronOccupancy = nullptr;

    G4DecayProducts* thePreAssignedDecayProducts = nullptr;

    // This void pointer is used by G4EventManager to maintain the
    // link between pre-assigned decay products and corresponding
    // primary particle
    G4PrimaryParticle* primaryParticle = nullptr;

    G4double theKineticEnergy = 0.0;

    mutable G4double theLogKineticEnergy = DBL_MAX;

    mutable G4double theBeta = -1.0;

    G4double theProperTime = 0.0;

    G4double theDynamicalMass = 0.0;

    G4double theDynamicalCharge = 0.0;

    G4double theDynamicalSpin = 0.0;

    G4double theDynamicalMagneticMoment = 0.0;

    G4double thePreAssignedDecayTime = -1.0;

    G4int verboseLevel = 1;

    G4int thePDGcode = 0;
};

#include "G4DynamicParticle.icc"

#endif
