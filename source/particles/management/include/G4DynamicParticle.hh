// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DynamicParticle.hh,v 1.2 1999-02-06 10:10:07 kurasige Exp $
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
//         modify thePreAssignedDecayProducts
//         add   void SetMomentum(G4ThreeVector &momentum)
//         add   void Set4Momentum(G4LorentzVector &momentum)
//         add   G4LorentzVector   Get4Momentum()
//         add   G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
//                                 G4LorentzVector &p4vector)
//      revised by Hisaya Kurashige, 19 Oct 1996
//         add    theKillProcess
//         add    theProperTime
//      revised by Hisaya Kurashige, 19 Feb 1997
//      revised by Hisaya Kurashige, 5  June 1998
//         remove    theKillProcess
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

class  G4VProcess;
class  G4DecayProducts;


class G4DynamicParticle 
{
  //  The dynamic particle is a class which contains the purely
  //  dynamic aspects of a moving particle. It also has a
  //  pointer to a G4ParticleDefinition object, which holds
  //  all the static information.

  public:
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

     ~G4DynamicParticle();

     G4DynamicParticle & operator=(const G4DynamicParticle &right);
     G4int operator==(const G4DynamicParticle &right) const;
     G4int operator!=(const G4DynamicParticle &right) const;

     inline void *operator new(size_t);
     inline void operator delete(void *aDynamicParticle);

     void DumpInfo() const;

     G4double GetMass() const;
     void     SetMass(G4double mass);
     // set/get dynamical mass
     // the dynamical mass is set to PDG mass in default

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

     G4double GetProperTime() const;
      //  Returns the current particle proper time
     void SetProperTime( G4double );
      //  Set the current particle Proper Time

     G4ParticleDefinition* GetDefinition() const;

     void SetDefinition(G4ParticleDefinition * aParticleDefinition);

     G4double GetTotalMomentum() const;
      //  Returns the module of the momentum vector

     G4double GetTotalEnergy() const;
      //  Returns the total energy of the particle

     const G4ThreeVector& GetPolarization() const;

     void SetPolarization(G4double polX, G4double polY, G4double polZ);

     G4double GetKineticEnergy() const;
      //  Returns the kinetic energy of a particle

     void SetKineticEnergy(G4double aEnergy);
      //  Sets the kinetic energy of a particle

     G4DecayProducts *GetPreAssignedDecayProducts() const;
     void SetPreAssignedDecayProducts(G4DecayProducts *aDecayProducts);

  private:
     G4double           theDynamicalMass;

     G4ThreeVector theMomentumDirection;
      //  The normalized momentum vector

     G4ParticleDefinition *theParticleDefinition;
      //  Contains the static information of this particle.

     G4ThreeVector thePolarization;

     G4double theKineticEnergy;

     G4double theProperTime;

     G4DecayProducts *thePreAssignedDecayProducts;

 private:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;


};

extern G4Allocator<G4DynamicParticle> aDynamicParticleAllocator;

// ------------------------
// Inlined operators
// ------------------------

inline void * G4DynamicParticle::operator new(size_t)
{
  void * aDynamicParticle;
  aDynamicParticle = (void *) aDynamicParticleAllocator.MallocSingle();
  return aDynamicParticle;
}

inline void G4DynamicParticle::operator delete(void * aDynamicParticle)
{
  aDynamicParticleAllocator.FreeSingle((G4DynamicParticle *) aDynamicParticle);
}

// ------------------------
// Inlined functions
// ------------------------
inline G4double G4DynamicParticle::GetMass() const
{
  return theDynamicalMass;
}

inline void G4DynamicParticle::SetMass(G4double newMass)
{
  theDynamicalMass = newMass;
}

inline const G4ThreeVector& G4DynamicParticle::GetMomentumDirection() const
{
  return theMomentumDirection;
}

inline G4ThreeVector G4DynamicParticle::GetMomentum() const
{
  G4double pModule = sqrt(theKineticEnergy*theKineticEnergy +
                        2*theKineticEnergy*theDynamicalMass);
  G4ThreeVector pMomentum(theMomentumDirection.x()*pModule,
                          theMomentumDirection.y()*pModule,
                          theMomentumDirection.z()*pModule);
  return pMomentum;
}

inline  G4LorentzVector  G4DynamicParticle::Get4Momentum() const
{       
  G4double mass      = theDynamicalMass;
  G4double energy    = theKineticEnergy;
  G4double momentum  = sqrt(energy*energy+2.0*mass*energy);
  G4LorentzVector    p4( theMomentumDirection.x()*momentum,
       			 theMomentumDirection.y()*momentum,
       			 theMomentumDirection.z()*momentum,
       			 energy+mass);
  return p4;
}

inline G4double G4DynamicParticle::GetTotalMomentum() const
{
  // The momentum is returned in energy equivalent.
  return sqrt((theKineticEnergy + 2.*theDynamicalMass)* theKineticEnergy);
}

inline G4ParticleDefinition* G4DynamicParticle::GetDefinition() const
{
  return theParticleDefinition;
}

inline const G4ThreeVector& G4DynamicParticle::GetPolarization() const
{
  return thePolarization;
}

inline G4double G4DynamicParticle::GetProperTime() const
{
  return theProperTime;
}

inline G4double G4DynamicParticle::GetTotalEnergy() const
{
  return (theKineticEnergy+theDynamicalMass);
}

inline G4double G4DynamicParticle::GetKineticEnergy() const
{
  return theKineticEnergy;
}

inline void G4DynamicParticle::SetMomentumDirection(const G4ThreeVector &aDirection)
{
  theMomentumDirection = aDirection;
}

inline void G4DynamicParticle::SetMomentumDirection(G4double px, G4double py, G4double pz)
{
  theMomentumDirection.setX(px);
  theMomentumDirection.setY(py);
  theMomentumDirection.setZ(pz);
}

inline void G4DynamicParticle::SetDefinition(G4ParticleDefinition * aParticleDefinition)
{
  theParticleDefinition = aParticleDefinition;
  theDynamicalMass = theParticleDefinition->GetPDGMass();
}

inline void G4DynamicParticle::SetPolarization(G4double polX, G4double polY, G4double polZ)
{
  thePolarization.setX(polX);
  thePolarization.setY(polY);
  thePolarization.setZ(polZ);
}

inline void G4DynamicParticle::SetKineticEnergy(G4double aEnergy)
{
  theKineticEnergy = aEnergy;
}

inline void G4DynamicParticle::SetProperTime(G4double atime)
{
  theProperTime = atime;
}

inline G4DecayProducts* G4DynamicParticle::GetPreAssignedDecayProducts() const
{ 
  return thePreAssignedDecayProducts;
}

inline void G4DynamicParticle::SetPreAssignedDecayProducts(G4DecayProducts* aDecayProducts)
{
 thePreAssignedDecayProducts = aDecayProducts;
}

inline 
void G4DynamicParticle::SetVerboseLevel(G4int value)
{
   verboseLevel = value;
}

inline 
G4int G4DynamicParticle::GetVerboseLevel() const
{
   return verboseLevel;
}

#endif










