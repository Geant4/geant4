// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEVector.hh,v 1.4 1999-12-15 16:42:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Gheisha friend class G4KinematicParticle  -- header file
// J.L. Chuma, TRIUMF, 22-Feb-1996
// last modified: H. Fesefeldt 18-November-1996

#ifndef G4HEVector_hh
#define G4HEVector_hh 1

#include "G4ParticleMomentum.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"

class G4HEVector
 {
  protected:
     G4double px;
     G4double py;
     G4double pz;
     G4double energy;
     G4double kineticEnergy;
     G4double mass;
     G4double charge;
     G4double timeOfFlight;
     G4int  side;
     G4bool flag;
     G4int  code;
     G4String particleName;
     G4String particleType;
     G4int    baryon;
     enum {NumberOfQuarkFlavor = 8};
     G4int theQuarkContent[NumberOfQuarkFlavor];
     G4int theAntiQuarkContent[NumberOfQuarkFlavor];

  public:

  G4HEVector(const G4DynamicParticle * aParticle);
  
  G4HEVector()
   {
     px               = 0.0;
     py               = 0.0;
     pz               = 0.0;
     energy           = 0.0;
     kineticEnergy    = 0.0;
     mass             = 0.0;
     charge           = 0.0;
     timeOfFlight     = 0.0;
     side             = 0;
     flag             = false;
     code             = 0;
     particleName     = "";
     particleType     = "";
     baryon           = 0;
   }


  G4HEVector( const G4HEVector & p )
   {
     px            = p.px;
     py            = p.py;
     pz            = p.pz;
     energy        = p.energy;
     kineticEnergy = p.kineticEnergy;
     mass          = p.mass;
     charge        = p.charge;
     timeOfFlight  = p.timeOfFlight;
     side          = p.side;
     flag          = p.flag;
     code          = p.code;
     particleName  = p.particleName;
     particleType  = p.particleType;
     baryon        = p.baryon;
   }


  G4HEVector & operator = ( const G4HEVector & p )
   {
     if ( this != &p ) {
       px            = p.px;
       py            = p.py;
       pz            = p.pz;
       energy        = p.energy;
       kineticEnergy = p.kineticEnergy;
       mass          = p.mass;
       charge        = p.charge;
       timeOfFlight  = p.timeOfFlight;
       side          = p.side;
       flag          = p.flag;
       code          = p.code;
       particleName  = p.particleName;
       particleType  = p.particleType;
       baryon        = p.baryon;
     }
     return *this;
   }

   ~G4HEVector(){ };

   G4String getParticleName(G4int code, G4int baryon);
 
   void setMomentum( const G4ParticleMomentum mom ); 

   void setMomentum( const G4ParticleMomentum * mom ); 

   void setMomentumAndUpdate( const G4ParticleMomentum mom );

   void setMomentumAndUpdate( const G4ParticleMomentum * mom );

   const G4ParticleMomentum getMomentum() const ;

   G4double getTotalMomentum();

   void setMomentum( G4double x, G4double y, G4double z);

   void setMomentumAndUpdate( G4double x, G4double y, G4double z );

   void setMomentum( G4double x, G4double y );

   void setMomentumAndUpdate( G4double x, G4double y );

   void setMomentum( G4double z );

   void setMomentumAndUpdate( G4double z );

   void setEnergy( G4double e ); 

   void setEnergyAndUpdate( G4double e );

   void setKineticEnergy( G4double ekin ); 

   void setKineticEnergyAndUpdate(G4double ekin); 

   G4double getEnergy(); 

   G4double getKineticEnergy(); 

   void setMass( G4double m ); 

   void setMassAndUpdate( G4double m );

   G4double getMass(); 

   void setCharge( G4double c ); 

   G4double getCharge(); 

   void setTOF( G4double t ); 

   G4double getTOF(); 

   void setSide( G4int s ); 

   G4int getSide(); 

   void setFlag( G4bool f ); 

   G4bool getFlag(); 

   void setCode( G4int c ); 

   G4int getCode(); 

   G4String getName();

   G4int getBaryonNumber();

   G4int getQuarkContent(G4int flavor);

   G4int getAntiQuarkContent(G4int flavor);

   void setZero();

   G4String getType();

   void Add( const G4HEVector & p1, const G4HEVector & p2 );

   void Sub( const G4HEVector & p1, const G4HEVector & p2 );

   void Lor( const G4HEVector & p1, const G4HEVector & p2 );

   G4double CosAng( const G4HEVector & p );

   G4double Ang(const G4HEVector & p );

   G4double Dot4( const G4HEVector & p1, const G4HEVector & p2);

   G4double Impu( const G4HEVector & p1, const G4HEVector & p2);

   void Add3( const G4HEVector & p1, const G4HEVector & p2);

   void Sub3( const G4HEVector & p1, const G4HEVector & p2);

   void Cross( const G4HEVector & p1, const G4HEVector & p2);

   G4double Dot( const G4HEVector & p1, const G4HEVector & p2);

   void Smul( const G4HEVector & p, G4double h);

   void SmulAndUpdate( const G4HEVector & p, G4double h);

   void Norz( const G4HEVector & p );

   G4double Length();

   void Exch( G4HEVector & p1);

   void Defs1( const G4HEVector & p1, const G4HEVector & p2);

   void Defs( const G4HEVector & p1, const G4HEVector & p2,
                    G4HEVector & my,       G4HEVector & mz );

   void Trac( const G4HEVector & p1, const G4HEVector & mx,
              const G4HEVector & my, const G4HEVector & mz);

   void setDefinition(G4String name);

   G4int FillQuarkContent();

   void Print( G4int L);
};

#endif


