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
//
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//      ---------------- G4DynamicParticle  ----------------
//      first implementation by Makoto Asai, 29 January 1996
//      revised by G.Cosmo, 29 February 1996
//      revised by Hisaya Kurashige, 24 July 1996
//      revised by Hisaya Kurashige, 19 Oct 1996
//      revised by Hisaya Kurashige, 19 Feb 1997
//                ------------------------
//      Add theDynamicCharge and theElectronOccupancy
//                             17 AUg. 1999   H.Kurashige  
//      Add thePreAssignedDecayTime   18 Jan. 2001 H.Kurashige
//      Added  MagneticMoment               Mar. 2007
//      Added GetLogKineticEnergy: the log-kinetic energy value is computed only 
 //     on demand if its stored value is not up-to-date. 15 March 2019 M. Novak
// ------------------------------------------------------------

#ifndef G4DynamicParticle_h
#define G4DynamicParticle_h 1

#include <cmath>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4Allocator.hh"
#include "G4LorentzVector.hh"
#include "G4Log.hh"

#include "G4ParticleMomentum.hh"
//  G4ParticleMomentum is "momentum direction" not "momentum vector"
//  The name is miss-leading so you should not use G4ParticleMomentum
//  and you are recommended to use G4ThreeVector instead

#include "G4ElectronOccupancy.hh"

class G4PrimaryParticle;
class G4DecayProducts;

class G4DynamicParticle 
{
  // Class Description
  //  The dynamic particle is a class which contains the purely
  //  dynamic aspects of a moving particle. It also has a
  //  pointer to a G4ParticleDefinition object, which holds
  //  all the static information.
  //

  public: // With Description
  //- constructors 
     G4DynamicParticle();

     G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aMomentumDirection,
                        G4double aKineticEnergy);
     G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aParticleMomentum);
     G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
                        const G4LorentzVector    &aParticleMomentum);
     G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
			G4double aTotalEnergy,
                        const G4ThreeVector &aParticleMomentum);
     G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
		       const G4ThreeVector& aMomentumDirection,
		       G4double aKineticEnergy,
		       const G4double dynamicalMass);

     G4DynamicParticle(const G4DynamicParticle &right);

  //- destructor
     ~G4DynamicParticle();

  //- operators
     G4DynamicParticle & operator=(const G4DynamicParticle &right);
     G4bool operator==(const G4DynamicParticle &right) const;
     G4bool operator!=(const G4DynamicParticle &right) const;

  //  new/delete operators are oberloded to use G4Allocator
     inline void *operator new(size_t);
     inline void operator delete(void *aDynamicParticle);

  //- Set/Get methods
 
     inline const G4ThreeVector& GetMomentumDirection() const;
      //  Returns the normalized direction of the momentum
     inline void SetMomentumDirection(const G4ThreeVector &aDirection);
      //  Sets the normalized direction of the momentum
     inline void SetMomentumDirection(G4double px, G4double py, G4double pz);
      //  Sets the normalized direction of the momentum by coordinates

     inline G4ThreeVector GetMomentum() const;
      //  Returns the current particle momentum vector
     void SetMomentum( const G4ThreeVector &momentum);
      //  set the current particle momentum vector

     inline G4LorentzVector Get4Momentum() const;
      //  Returns the current particle energy-momentum 4vector
     void Set4Momentum( const G4LorentzVector &momentum);
      //  Set the current particle energy-momentum 4vector

     inline G4double GetTotalMomentum() const;
      //  Returns the module of the momentum vector
     inline G4double GetTotalEnergy() const;
      //  Returns the total energy of the particle

     inline G4double GetKineticEnergy() const;
      //  Returns the kinetic energy of a particle
     inline G4double GetLogKineticEnergy() const;
      //  Returns:
      //   - natural logarithm of the particle kinetic energy (E_k) if E_k > 0
      //   - LOG_EKIN_MIN otherwise
     inline void SetKineticEnergy(G4double aEnergy);
      //  Sets the kinetic energy of a particle

     inline G4double GetProperTime() const;
      //  Returns the current particle proper time
     inline void SetProperTime( G4double );
      //  Set the current particle Proper Time

     inline const G4ThreeVector& GetPolarization() const;
     inline void SetPolarization(const G4ThreeVector&);
     inline void SetPolarization(G4double polX, G4double polY, G4double polZ);
      //   Set/Get polarization vector       

     inline G4double GetMass() const;
     inline void     SetMass(G4double mass);
     // set/get dynamical mass
     // the dynamical mass is set to PDG mass in default

     inline G4double GetCharge() const;
     inline void     SetCharge(G4double charge);
     inline void     SetCharge(G4int    chargeInUnitOfEplus);
     // set/get dynamical charge 
     // the dynamical mass is set to PDG charge in default

     inline G4double GetSpin() const;
     inline void     SetSpin(G4double spin);
     inline void     SetSpin(G4int    spinInUnitOfHalfInteger);
     // set/get dynamical spin
     // the dynamical spin is set to PDG spin in default

     inline G4double GetMagneticMoment() const;
     inline void     SetMagneticMoment(G4double magneticMoment);
     // set/get dynamical MagneticMoment  
     // the dynamical mass is set to PDG MagneticMoment in default

     inline const G4ElectronOccupancy* GetElectronOccupancy() const;
     // Get electron occupancy 
     // ElectronOccupancy is valid only if the particle is ion
     inline G4int  GetTotalOccupancy() const;
     inline G4int  GetOccupancy(G4int orbit) const;
     inline void   AddElectron(G4int orbit, G4int number = 1);
     inline void   RemoveElectron(G4int orbit, G4int number = 1);
  
     inline const G4ParticleDefinition* GetParticleDefinition() const;
     void SetDefinition(const G4ParticleDefinition * aParticleDefinition);
     //   Set/Get particle definition  
     //  following method of GetDefinition remains 
     //  because of backward compatiblity. It will be removed in future 
     inline G4ParticleDefinition* GetDefinition() const;
     
     inline const G4DecayProducts *GetPreAssignedDecayProducts() const;
     inline void SetPreAssignedDecayProducts(G4DecayProducts *aDecayProducts);
      //   Set/Get pre-assigned decay channel

     inline G4double GetPreAssignedDecayProperTime() const;
     inline void SetPreAssignedDecayProperTime(G4double);
      //   Set/Get pre-assigned proper time when the particle will decay
   
     //- print out information
     void DumpInfo(G4int mode= 0) const;
     //    mode 0 : default )(minimum)
     //    mode 1 : 0 + electron occupancy

 protected:
     void      AllocateElectronOccupancy(); 
     G4double  GetElectronMass() const;

 public:  // With Description
   inline void  SetVerboseLevel(G4int value);
   inline G4int GetVerboseLevel() const;
   // Set/Get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

   inline void SetPrimaryParticle(G4PrimaryParticle* p);
   inline void SetPDGcode(G4int c);

 public: // With Description
   inline G4PrimaryParticle* GetPrimaryParticle() const;
   // Return the pointer to the corresponding G4PrimaryParticle object
   // if this particle is a primary particle OR is defined as a pre-assigned
   // decay product. Otherwise return null.

   inline G4int GetPDGcode() const;
   // Return the PDG code of this particle. If the particle is known to Geant4
   // its PDG code defined in G4ParticleDefinition is returned. If it is unknown
   // (i.e. PDG code in G4ParticleDefinition is 0), PDG code defined in the
   // corresponding primary particle or pre-assigned decay product will be
   // returned if available. Otherwise (e.g. for geantino) returns 0.

 private:

   G4ThreeVector theMomentumDirection;
   //  The normalized momentum vector

   G4ThreeVector thePolarization;

   const G4ParticleDefinition *theParticleDefinition;
   //  Contains the static information of this particle.

   G4ElectronOccupancy* theElectronOccupancy;          
  
   G4DecayProducts *thePreAssignedDecayProducts;

   G4PrimaryParticle* primaryParticle;
   // This void pointer is used by G4EventManager to maintain the
   // link between pre-assigned decay products and corresponding
   // primary particle.

   G4double theKineticEnergy;

   mutable G4double theLogKineticEnergy;

   G4double theProperTime;

   G4double theDynamicalMass;

   G4double theDynamicalCharge;

   G4double theDynamicalSpin;

   G4double theDynamicalMagneticMoment;

   G4double thePreAssignedDecayTime;

   G4int verboseLevel;

   G4int thePDGcode;
};

#include "G4DynamicParticle.icc"

#endif
