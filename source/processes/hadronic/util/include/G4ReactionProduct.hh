// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReactionProduct.hh,v 1.2 1999-12-15 14:53:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // J.L. Chuma, TRIUMF, 31-Oct-1996
 // last modified: 19-Dec-1996
 // modified by J.L.Chuma, 24-Jul-1997   to include total momentum
 // inluded operator *, and some minor modifications.
 // modified by H.P.Wellisch to add functionality needed by string models,
 // cascade and Nucleus. (Mon Mar 16 1998) 
 
#ifndef G4ReactionProduct_h
#define G4ReactionProduct_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
 
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
       result.SetTotalEnergy(sqrt(result.GetMass()*result.GetMass()+
                                  result.GetMomentum()*result.GetMomentum()));
       return result;
     }

 public:
    G4ReactionProduct();
    
    G4ReactionProduct( G4ParticleDefinition *aParticleDefinition );
    
    ~G4ReactionProduct() {}
    
    G4ReactionProduct( const G4ReactionProduct &right );
    
    G4ReactionProduct &operator= ( const G4ReactionProduct &right );
    
    G4ReactionProduct &operator= ( const G4DynamicParticle &right );
    
    inline G4bool operator== ( const G4ReactionProduct &right ) const
    { return ( this == (G4ReactionProduct*) &right ); }
    
    inline G4bool operator!= ( const G4ReactionProduct &right ) const
    { return ( this != (G4ReactionProduct*) &right ); }
    
    inline G4ParticleDefinition *GetDefinition() const
    { return theParticleDefinition; }
    
    void SetDefinition( G4ParticleDefinition *aParticleDefinition );
   
    void SetDefinitionAndUpdateE( G4ParticleDefinition *aParticleDefinition );
      
    void SetMomentum( const G4double x, const G4double y, const G4double z );
    
    void SetMomentum( const G4double x, const G4double y );
    
    void SetMomentum( const G4double z );

    inline void SetMomentum( const G4ThreeVector &m )
    { momentum = m; }
    
    inline G4ThreeVector GetMomentum() const
    { return momentum; }
    
    inline G4double GetTotalMomentum() const
    { return sqrt(abs(kineticEnergy*(totalEnergy+mass))); }
    
    inline G4double GetTotalEnergy() const
    { return totalEnergy; }
    
    inline void SetKineticEnergy( const G4double e )
    {
      kineticEnergy = e;
      totalEnergy = kineticEnergy + mass;
    }
    
    inline G4double GetKineticEnergy() const
    { return kineticEnergy; }

    inline void SetTotalEnergy( const G4double e )
    {
      totalEnergy = e;
      kineticEnergy = totalEnergy - mass;
    }
    
    inline void SetMass( const G4double m )
    { mass = m; }
    
    inline G4double GetMass() const
    { return mass; }
    
    inline void SetTOF( const G4double t )
    { timeOfFlight = t; }
    
    inline G4double GetTOF() const
    { return timeOfFlight; }
    
    inline void SetSide( const G4int s )
    { side = s; }
    
    inline G4int GetSide() const
    { return side; }
    
    inline void SetNewlyAdded( const G4bool f )
    { NewlyAdded = f; }
    
    inline G4bool GetNewlyAdded() const
    { return NewlyAdded; }
    
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
    
    inline G4ThreeVector GetPositionInNucleus() const {return positionInNucleus; }
    inline G4double GetXPositionInNucleus() const { return positionInNucleus.x(); }
    inline G4double GetYPositionInNucleus() const { return positionInNucleus.y(); }
    inline G4double GetZPositionInNucleus() const { return positionInNucleus.z(); }
    
    inline void SetFormationTime(G4double aTime) { formationTime = aTime; }
    
    inline G4double GetFormationTime() const { return formationTime; }
    
    inline void HasInitialStateParton(G4bool aFlag) { hasInitialStateParton = aFlag; }
    
    inline G4bool HasInitialStateParton() const { return hasInitialStateParton; }
    
 private:
    
    G4ParticleDefinition *theParticleDefinition;
    
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

    // NewlyAdded refers to particles added by "nuclear excitation", or as
    //  "black track" particles, or as deuterons, tritons, and alphas
    G4bool NewlyAdded;
 };
 
#endif
 
