// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DynamicParticle.hh,v 1.9 2001-01-18 12:46:05 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
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
// ------------------------------------------------------------

#ifndef G4DynamicParticle_h
#define G4DynamicParticle_h 1


#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4Allocator.hh"
#include "G4LorentzVector.hh"

#include "G4ParticleMomentum.hh"
//  G4ParticleMomentum is "momentum direction" not "momentum vector"
//  The name is miss-leading so you should not use G4ParticleMomentum
//  and you are recommended to use G4ThreeVector instead

#include "G4ElectronOccupancy.hh"

class  G4VProcess;
class  G4DecayProducts;


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

     G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aMomentumDirection,
                        G4double aKineticEnergy);
     G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aParticleMomentum);
     G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                        const G4LorentzVector    &aParticleMomentum);
     G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
			G4double aTotalEnergy,
                        const G4ThreeVector &aParticleMomentum);

     G4DynamicParticle(const G4DynamicParticle &right);

  //- destructor
     ~G4DynamicParticle();

  //- operators
     G4DynamicParticle & operator=(const G4DynamicParticle &right);
     G4int operator==(const G4DynamicParticle &right) const;
     G4int operator!=(const G4DynamicParticle &right) const;

  //  new/delete operators are oberloded to use G4Allocator
     inline void *operator new(size_t);
     inline void operator delete(void *aDynamicParticle);

  //- Set/Get methods
 
     const G4ThreeVector& GetMomentumDirection() const;
      //  Returns the normalized direction of the momentum
     void SetMomentumDirection(const G4ThreeVector &aDirection);
      //  Sets the normalized direction of the momentum
     void SetMomentumDirection(G4double px, G4double py, G4double pz);
      //  Sets the normalized direction of the momentum by coordinates

     G4ThreeVector GetMomentum() const;
      //  Returns the current particle momentum vector
     void SetMomentum( const G4ThreeVector &momentum);
      //  set the current particle momentum vector

     G4LorentzVector Get4Momentum() const;
      //  Returns the current particle energy-momentum 4vector
     void Set4Momentum( const G4LorentzVector &momentum);
      //  Set the current particle energy-momentum 4vector


     G4double GetTotalMomentum() const;
      //  Returns the module of the momentum vector
     G4double GetTotalEnergy() const;
      //  Returns the total energy of the particle

     G4double GetKineticEnergy() const;
      //  Returns the kinetic energy of a particle
     void SetKineticEnergy(G4double aEnergy);
      //  Sets the kinetic energy of a particle


     G4double GetProperTime() const;
      //  Returns the current particle proper time
     void SetProperTime( G4double );
      //  Set the current particle Proper Time


     const G4ThreeVector& GetPolarization() const;
     void SetPolarization(G4double polX, G4double polY, G4double polZ);
      //   Set/Get polarization vector       


     G4double GetMass() const;
     void     SetMass(G4double mass);
     // set/get dynamical mass
     // the dynamical mass is set to PDG mass in default


     G4double GetCharge() const;
     void     SetCharge(G4double charge);
     void     SetCharge(G4int    chargeInUnitOfEplus);
     // set/get dynamical charge 
     // the dynamical mass is set to PDG charge in default

     const G4ElectronOccupancy* GetElectronOccupancy() const;
     // Get electron occupancy 
     // ElectronOccupancy is valid only if the particle is ion
     G4int  GetTotalOccupancy() const;
     G4int  GetOccupancy(G4int orbit) const;
     void   AddElectron(G4int orbit, G4int number = 1);
     void   RemoveElectron(G4int orbit, G4int number = 1);
 
     G4ParticleDefinition* GetDefinition() const;
     void SetDefinition(G4ParticleDefinition * aParticleDefinition);
      //   Set/Get particle definition  

     
     const G4DecayProducts *GetPreAssignedDecayProducts() const;
     void SetPreAssignedDecayProducts(G4DecayProducts *aDecayProducts);
      //   Set/Get pre-assigned decay channel

     G4double GetPreAssignedDecayProperTime() const;
     void SetPreAssignedDecayProperTime(G4double);
      //   Set/Get pre-assigned proper time when the particle will decay
 
   
  //- print out information
     void DumpInfo(G4int mode= 0) const;
     //    mode 0 : default )(minimum)
     //    mode 1 : 0 + electron occupancy

   private:
     void      AllocateElectronOccupancy(); 
     G4double  GetElectronMass() const;

   private:
     G4ThreeVector theMomentumDirection;
      //  The normalized momentum vector

     G4ParticleDefinition *theParticleDefinition;
      //  Contains the static information of this particle.

     G4ThreeVector thePolarization;

     G4double theKineticEnergy;

     G4double theProperTime;

     G4double theDynamicalMass;

     G4double theDynamicalCharge;

     G4ElectronOccupancy* theElectronOccupancy;          
  
     G4DecayProducts *thePreAssignedDecayProducts;

     G4double thePreAssignedDecayTime;

 private:
   G4int verboseLevel;
 
 public:  // With Description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // Set/Get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More


};

#include "G4DynamicParticle.icc"

#endif










