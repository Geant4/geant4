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
// J.L. Chuma, TRIUMF, 31-Oct-1996
// last modified: 19-Dec-1996
// modified by J.L.Chuma, 24-Jul-1997   to include total momentum
// inluded operator *, and some minor modifications.
// modified by H.P.Wellisch to add functionality needed by string models,
// cascade and Nucleus. (Mon Mar 16 1998) 
// M. Kelsey 29-Aug-2011 -- Use G4Allocator model to avoid memory churn.
 
#ifndef G4ReactionProduct_h
#define G4ReactionProduct_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4DynamicParticle.hh"
#include "G4HadProjectile.hh"
#include "G4HadronicException.hh"

class G4ReactionProduct;

// To support better memory management and reduced fragmentation
//
#if defined G4HADRONIC_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4ReactionProduct>*& aRPAllocator();
#else
  extern G4DLLIMPORT G4Allocator<G4ReactionProduct>*& aRPAllocator();
#endif

class G4ReactionProduct
{
    friend G4ReactionProduct operator+(
     const G4ReactionProduct & p1, const G4ReactionProduct &p2 );
    
    friend G4ReactionProduct operator-(
     const G4ReactionProduct & p1, const G4ReactionProduct &p2 );

    friend G4ReactionProduct operator*(
     const G4double aDouble, const G4ReactionProduct &p2 )
     {
       G4ReactionProduct result;
       result.SetMomentum(aDouble*p2.GetMomentum());
       result.SetMass(p2.GetMass());
       result.SetTotalEnergy(std::sqrt(result.GetMass()*result.GetMass()+
                                  result.GetMomentum()*result.GetMomentum()));
       return result;
     }

 public:
    G4ReactionProduct();
    
    G4ReactionProduct(const G4ParticleDefinition *aParticleDefinition );

    ~G4ReactionProduct() {}
    
    G4ReactionProduct( const G4ReactionProduct &right );

    // Override new and delete for use with G4Allocator
    inline void* operator new(size_t) {
      if (!aRPAllocator()) aRPAllocator() = new G4Allocator<G4ReactionProduct>  ;
      return (void *)aRPAllocator()->MallocSingle();
    }
#ifdef __IBMCPP__
    inline void* operator new(size_t, void *p) {
      return p;
    }
#endif
    inline void operator delete(void* aReactionProduct) {
      aRPAllocator()->FreeSingle((G4ReactionProduct*)aReactionProduct);
    }

    G4ReactionProduct &operator= ( const G4ReactionProduct &right );
    
    G4ReactionProduct &operator= ( const G4DynamicParticle &right );
    
    G4ReactionProduct &operator= ( const G4HadProjectile &right );

    inline G4bool operator== ( const G4ReactionProduct &right ) const
    { return ( this == (G4ReactionProduct*) &right ); }
    
    inline G4bool operator!= ( const G4ReactionProduct &right ) const
    { return ( this != (G4ReactionProduct*) &right ); }
    
    inline const G4ParticleDefinition* GetDefinition() const
    { return theParticleDefinition; }

    void SetDefinition(const G4ParticleDefinition* aParticleDefinition );
   
    void SetDefinitionAndUpdateE(const G4ParticleDefinition* aParticleDefinition );
      
    void SetMomentum( const G4double x, const G4double y, const G4double z );
    
    void SetMomentum( const G4double x, const G4double y );
    
    void SetMomentum( const G4double z );

    inline void SetMomentum( const G4ThreeVector &mom )
    { momentum = mom; }
    
    inline G4ThreeVector GetMomentum() const
    { return momentum; }
    
    inline G4double GetTotalMomentum() const
    { return std::sqrt(std::abs(kineticEnergy*(totalEnergy+mass))); }
    
    inline G4double GetTotalEnergy() const
    { return totalEnergy; }
    
    inline void SetKineticEnergy( const G4double en )
    {
      kineticEnergy = en;
      totalEnergy = kineticEnergy + mass;
    }
    
    inline G4double GetKineticEnergy() const
    { return kineticEnergy; }

    inline void SetTotalEnergy( const G4double en )
    {
      totalEnergy = en;
      kineticEnergy = totalEnergy - mass;
    }
    
    inline void SetMass( const G4double mas )
    { mass = mas; }
    
    inline G4double GetMass() const
    { return mass; }
    
    inline void SetTOF( const G4double t )
    { timeOfFlight = t; }
    
    inline G4double GetTOF() const
    { return timeOfFlight; }
    
    inline void SetSide( const G4int sid )
    { side = sid; }
    
    inline G4int GetSide() const
    { return side; }
    
    inline void SetCreatorModelID( const G4int mod )
    { theCreatorModel = mod; }
    
    inline G4int GetCreatorModelID() const
    { return theCreatorModel; }

    inline const G4ParticleDefinition* GetParentResonanceDef() const
    { return theParentResonanceDef; }

    inline void SetParentResonanceDef( const G4ParticleDefinition* parentDef )
    { theParentResonanceDef = parentDef; }

    inline G4int GetParentResonanceID() const { return theParentResonanceID; }

    inline void SetParentResonanceID ( const G4int parentID )
    { theParentResonanceID = parentID; }
    
    inline void SetNewlyAdded( const G4bool f )
    { NewlyAdded = f; }
    
    inline G4bool GetNewlyAdded() const
    { return NewlyAdded; }
    
    inline void SetMayBeKilled( const G4bool f )
    { MayBeKilled = f; }
    
    inline G4bool GetMayBeKilled() const
    { return MayBeKilled; }

    void SetZero();
    
    void Lorentz( const G4ReactionProduct &p1, const G4ReactionProduct &p2 );
    
    G4double Angle( const G4ReactionProduct &p ) const;
    
    inline void SetPositionInNucleus(G4double x, G4double y, G4double z)
     {
       positionInNucleus.setX(x);
       positionInNucleus.setY(y);
       positionInNucleus.setZ(z);
     }
    
    inline void SetPositionInNucleus( G4ThreeVector & aPosition )
     {
       positionInNucleus = aPosition;
     }
    
    inline G4ThreeVector GetPositionInNucleus() const { return positionInNucleus; }
    inline G4double GetXPositionInNucleus() const { return positionInNucleus.x(); }
    inline G4double GetYPositionInNucleus() const { return positionInNucleus.y(); }
    inline G4double GetZPositionInNucleus() const { return positionInNucleus.z(); }
    
    inline void SetFormationTime(G4double aTime) { formationTime = aTime; }
    
    inline G4double GetFormationTime() const { return formationTime; }
    
    inline void HasInitialStateParton(G4bool aFlag) { hasInitialStateParton = aFlag; }
    
    inline G4bool HasInitialStateParton() const { return hasInitialStateParton; }
 
private:
    
    const G4ParticleDefinition *theParticleDefinition;
    
    // for use with string models and cascade.
    G4ThreeVector positionInNucleus;
    G4double formationTime;
    G4bool hasInitialStateParton;
    
    // mass is included here, since pseudo-particles are created with masses different
    // than the standard particle masses, and we are not allowed to create particles
    G4double mass;
    
    G4ThreeVector momentum;
    
    G4double totalEnergy;
    G4double kineticEnergy;
    
    G4double timeOfFlight;
    
    //  side refers to how the particles are distributed in the
    //  forward (+) and backward (-) hemispheres in the center of mass system
    G4int side;

    G4int theCreatorModel;

    const G4ParticleDefinition* theParentResonanceDef = nullptr;
    G4int theParentResonanceID;

    // NewlyAdded refers to particles added by "nuclear excitation", or as
    //  "black track" particles, or as deuterons, tritons, and alphas
    G4bool NewlyAdded;
    G4bool MayBeKilled;
};
 
#endif
 
