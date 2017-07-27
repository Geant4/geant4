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
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4GHEKinematicsVector utility class ------
//                   by Larry Felawka (TRIUMF), March 1997
//                     E-mail: felawka@alph04.triumf.ca
// ************************************************************
//-----------------------------------------------------------------------------

// Store, Retrieve and manipulate particle data.
// Based on "G4GHEVector" class of H. Fesefeldt.

#ifndef G4GHEKinematicsVector_h
#define G4GHEKinematicsVector_h 1

#include "G4ios.hh"

#include <CLHEP/Units/PhysicalConstants.h>

class G4GHEKinematicsVector
 {
 public:
  inline
   G4GHEKinematicsVector()
   {
     momentum.setX(  0.0 );
     momentum.setY(  0.0 );
     momentum.setZ(  0.0 );
     energy        = 0.0;
     kineticEnergy = 0.0;
     mass          = 0.0;
     charge        = 0.0;
     timeOfFlight  = 0.0;
     side          = 0;
     flag          = false;
     code          = 0;
     particleDef   = NULL;
   }

  ~G4GHEKinematicsVector() {}

  inline
   G4GHEKinematicsVector( const G4GHEKinematicsVector & p )
   {
     momentum.setX( p.momentum.x() );
     momentum.setY( p.momentum.y() );
     momentum.setZ( p.momentum.z() );
     energy        = p.energy;
     kineticEnergy = p.kineticEnergy;
     mass          = p.mass;
     charge        = p.charge;
     timeOfFlight  = p.timeOfFlight;
     side          = p.side;
     flag          = p.flag;
     code          = p.code;
     particleDef   = p.particleDef;
   }

  inline
   G4GHEKinematicsVector & operator = ( const G4GHEKinematicsVector & p )
   {
     if (this != &p)
     {
        momentum.setX( p.momentum.x() );
        momentum.setY( p.momentum.y() );
        momentum.setZ( p.momentum.z() );
        energy        = p.energy;
        kineticEnergy = p.kineticEnergy;
        mass          = p.mass;
        charge        = p.charge;
        timeOfFlight  = p.timeOfFlight;
        side          = p.side;
        flag          = p.flag;
        code          = p.code;
        particleDef   = p.particleDef;
     }
    return *this;
   }

  inline
   void SetMomentum( G4ParticleMomentum mom ) { momentum = mom; return; };

  inline
   void SetMomentumAndUpdate( G4ParticleMomentum mom )
   {
     momentum      = mom;
     energy        = std::sqrt(mass*mass + momentum.mag2());
     kineticEnergy = std::max(0.,energy - mass);
     return;
   }

  inline const
  G4ParticleMomentum GetMomentum() const { return momentum; }

  inline
   void SetMomentum( G4double x, G4double y, G4double z)
   { 
     momentum.setX( x );
     momentum.setY( y );
     momentum.setZ( z );
     return;
   } 

  inline
   void SetMomentumAndUpdate( G4double x, G4double y, G4double z )
   {
     momentum.setX( x );
     momentum.setY( y );
     momentum.setZ( z );
     energy        = std::sqrt(mass*mass + momentum.mag2());
     kineticEnergy = std::max(0.,energy-mass);
     return;
   }

  inline
   void SetMomentum( G4double x, G4double y )
   {
     momentum.setX( x );
     momentum.setY( y );
     return;
   }

  inline
   void SetMomentumAndUpdate( G4double x, G4double y )
   {
     momentum.setX( x );
     momentum.setY( y );
     energy = std::sqrt(mass*mass + momentum.mag2());
     kineticEnergy = std::max(0.,energy-mass);
     return;
   }

  inline
   void SetMomentum( G4double z )
   {
     momentum.setZ( z );
     return;
   }

  inline
   void SetMomentumAndUpdate( G4double z )
   {
     momentum.setZ( z );
     energy = std::sqrt(mass*mass + momentum.mag2());
     kineticEnergy = std::max(0.,energy-mass);
     return;
   }

  inline 
   void SetEnergy( G4double e ) { energy = e; return; }

  inline
   void SetEnergyAndUpdate( G4double e )
   {
     if (e <= mass)
       { 
         energy        = mass;
         kineticEnergy = 0.;
         momentum.setX(  0.);
         momentum.setY(  0.);
         momentum.setZ(  0.);
       }
     else
       {
         energy = e;
         kineticEnergy   = energy - mass;
         G4double momold = momentum.mag();
         G4double momnew = std::sqrt(energy*energy - mass*mass);
         if (momold == 0.)
           {
             G4double cost = 1.0- 2.0*G4UniformRand();
             G4double sint = std::sqrt(1. - cost*cost);
             G4double phi  = CLHEP::twopi* G4UniformRand();
             momentum.setX( momnew * sint * std::cos(phi));
             momentum.setY( momnew * sint * std::sin(phi));
             momentum.setZ( momnew * cost);
           }
         else
           {
             momnew /= momold;
             momentum.setX(momentum.x()*momnew);
             momentum.setY(momentum.y()*momnew);
             momentum.setZ(momentum.z()*momnew);
           }
       }    
     return;
   }

  inline 
   void SetKineticEnergy( G4double ekin ) { kineticEnergy = ekin; return; }

  inline 
   void SetKineticEnergyAndUpdate(G4double ekin) 
   {
     if (ekin <= 0.)
       { 
         energy        = mass;
         kineticEnergy = 0.;
         momentum.setX(  0.);
         momentum.setY(  0.);
         momentum.setZ(  0.);
       }
     else
       {
         energy = ekin + mass;
         kineticEnergy   = ekin;
         G4double momold = momentum.mag();
         G4double momnew = std::sqrt(energy*energy - mass*mass);
         if (momold == 0.)
           {
             G4double cost = 1.0-2.0*G4UniformRand();
             G4double sint = std::sqrt(1. - cost*cost);
             G4double phi  = CLHEP::twopi* G4UniformRand();
             momentum.setX( momnew * sint * std::cos(phi));
             momentum.setY( momnew * sint * std::sin(phi));
             momentum.setZ( momnew * cost);
           }
         else
           {
             momnew /= momold;
             momentum.setX(momentum.x()*momnew);
             momentum.setY(momentum.y()*momnew);
             momentum.setZ(momentum.z()*momnew);
           }
       }    
     return;
   }

  inline
  G4double GetEnergy() {return energy;}

  inline
  G4double GetKineticEnergy() {return kineticEnergy;}

  inline
  void SetMass( G4double mas ) { mass = mas; return; }

  inline
  void SetMassAndUpdate( G4double mas )
  {
    kineticEnergy = std::max(0., energy - mas);
    mass = mas;
    energy = kineticEnergy + mass;
    G4double momnew = std::sqrt(std::max(0., energy*energy - mass*mass));
    if ( momnew == 0.0) 
       {
         momentum.setX( 0.0 );
         momentum.setY( 0.0 );
         momentum.setZ( 0.0 );
       }
    else
       {
         G4double momold = momentum.mag();
         if (momold == 0.)
            { 
              G4double cost = 1.-2.*G4UniformRand();
              G4double sint = std::sqrt(1.-cost*cost);
              G4double phi  = CLHEP::twopi*G4UniformRand();
              momentum.setX( momnew*sint*std::cos(phi));
              momentum.setY( momnew*sint*std::sin(phi));
              momentum.setZ( momnew*cost);
            }
         else
            {
              momnew /= momold;
              momentum.setX( momentum.x()*momnew );
              momentum.setY( momentum.y()*momnew );
              momentum.setZ( momentum.z()*momnew );
            }
       }     
    return;
  }    

  inline
  G4double GetMass() { return mass; }

  inline
  void SetCharge( G4double c ) { charge = c; return; }

  inline
  G4double GetCharge() {return charge; }

  inline
  void SetTOF( G4double t ) { timeOfFlight = t; return; }

  inline
  G4double GetTOF() { return timeOfFlight; }

  inline
  void SetSide( G4int sid ) { side = sid; return; }

  inline
  G4int GetSide() { return side; }

  inline
  void setFlag( G4bool f ) { flag = f; return; }

  inline
  G4bool getFlag() { return flag; }

  inline
  void SetCode( G4int c ) { code = c; return; }

  inline
  void SetParticleDef( G4ParticleDefinition * c ) { particleDef = c; return; }

  inline
  G4int GetCode() { return code; } 

  inline
  G4ParticleDefinition * GetParticleDef() { return particleDef; } 

  inline
  void SetZero()
   {
     momentum.setX(  0.0 );
     momentum.setY(  0.0 );
     momentum.setZ(  0.0 );
     energy        = 0.0;
     kineticEnergy = 0.0;
     mass          = 0.0;
     charge        = 0.0;
     timeOfFlight  = 0.0;
     side          = 0;
     flag          = false;
     code          = 0;
     particleDef   = NULL;
   }

  inline
   void Add( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2 )
   {
     momentum = p1.momentum + p2.momentum;
     energy = p1.energy + p2.energy;
     G4double b = energy*energy - momentum.mag2();
     if( b < 0 )
       mass = -1. * std::sqrt( -b );
     else
       mass = std::sqrt( b );
     kineticEnergy = std::max(0.,energy - mass);
     charge        = p1.charge + p2.charge;
     code          = p1.code   + p2.code;
     particleDef   = p1.particleDef;
   }

  inline
   void Sub( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2 )
   {
     momentum = p1.momentum - p2.momentum;
     energy = p1.energy - p2.energy;
     G4double b = energy*energy - momentum.mag2();
     if( b < 0 )
       mass = -1. * std::sqrt( -b );
     else
       mass = std::sqrt( b );
     kineticEnergy = std::max(0.,energy - mass);
     charge        = p1.charge - p2.charge;
     code          = p1.code   - p2.code;
     particleDef   = p1.particleDef;
   }

  inline
   void Lor( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2 )
   {
     G4double a;
     a = ( p1.momentum.dot(p2.momentum)/(p2.energy+p2.mass) - p1.energy ) / p2.mass;
     momentum.setX( p1.momentum.x()+a*p2.momentum.x() );
     momentum.setY( p1.momentum.y()+a*p2.momentum.y() );
     momentum.setZ( p1.momentum.z()+a*p2.momentum.z() );
     energy = std::sqrt( sqr(p1.mass) + momentum.mag2() );
     mass = p1.mass;
     kineticEnergy = std::max(0.,energy - mass);
     timeOfFlight  = p1.timeOfFlight;
     side          = p1.side;
     flag          = p1.flag;
     code          = p1.code;
     particleDef   = p1.particleDef;
   }

  inline
   G4double CosAng( const G4GHEKinematicsVector & p )
   {
     G4double a = std::sqrt( momentum.mag2() * p.momentum.mag2() );
     if( a != 0.0 ) 
       {
         a = (momentum.x()*p.momentum.x() +
              momentum.y()*p.momentum.y() +
              momentum.z()*p.momentum.z()) / a;
         if( std::fabs(a) > 1.0 ) a<0.0 ? a=-1.0 : a=1.0;
       }
     return a;
   }
  inline
   G4double Ang(const G4GHEKinematicsVector & p )
   {
     G4double a = std::sqrt( momentum.mag2() * p.momentum.mag2() );
     if( a != 0.0 ) 
       {
         a = (momentum.x()*p.momentum.x() +
              momentum.y()*p.momentum.y() +
              momentum.z()*p.momentum.z()) / a;
         if( std::fabs(a) > 1.0 ) a<0.0 ? a=-1.0 : a=1.0;
       }
     return std::acos(a);
   }      

  inline
   G4double Dot4( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     return (   p1.energy       * p2.energy
              - p1.momentum.x() * p2.momentum.x()
              - p1.momentum.y() * p2.momentum.y()
              - p1.momentum.z() * p2.momentum.z() );
   } 
  
  inline 
   G4double Impu( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     return ( - sqr( p1.energy      - p2.energy)  
              + sqr(p1.momentum.x() - p2.momentum.x())
              + sqr(p1.momentum.y() - p2.momentum.y()) 
              + sqr(p1.momentum.z() - p2.momentum.z()) );
   }

  inline 
   void Add3( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     momentum.setX( p1.momentum.x() + p2.momentum.x());
     momentum.setY( p1.momentum.y() + p2.momentum.y());
     momentum.setZ( p1.momentum.z() + p2.momentum.z());    
     return;
   }

  inline
   void Sub3( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     momentum.setX( p1.momentum.x() - p2.momentum.x());
     momentum.setY( p1.momentum.y() - p2.momentum.y());
     momentum.setZ( p1.momentum.z() - p2.momentum.z());
     return;
   }

  inline
   void Cross( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     G4double px, py, pz;
     px = p1.momentum.y() * p2.momentum.z() - p1.momentum.z() * p2.momentum.y();
     py = p1.momentum.z() * p2.momentum.x() - p1.momentum.x() * p2.momentum.z();
     pz = p1.momentum.x() * p2.momentum.y() - p1.momentum.y() * p2.momentum.x();
     momentum.setX( px );
     momentum.setY( py );
     momentum.setZ( pz ); 
     return;
   }  

  inline 
   G4double Dot( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     return (   p1.momentum.x() * p2.momentum.x()
              + p1.momentum.y() * p2.momentum.y()
              + p1.momentum.z() * p2.momentum.z() );
   }

  inline
   void Smul( const G4GHEKinematicsVector & p, G4double h)
   {
     momentum.setX( h * p.momentum.x());
     momentum.setY( h * p.momentum.y());
     momentum.setZ( h * p.momentum.z());
     return;
   }   

  inline
   void SmulAndUpdate( const G4GHEKinematicsVector & p, G4double h)
   {
     momentum.setX( h * p.momentum.x());
     momentum.setY( h * p.momentum.y());
     momentum.setZ( h * p.momentum.z());
     mass          = p.mass;
     energy        = std::sqrt(momentum.mag2() + mass*mass);
     kineticEnergy = energy - mass;
     charge        = p.charge;
     timeOfFlight  = p.timeOfFlight;
     side          = p.side;
     flag          = p.flag;
     code          = p.code;
     particleDef   = p.particleDef;
     return;
   }
 
  inline
   void Norz( const G4GHEKinematicsVector & p )
   {
     G4double a =   p.momentum.mag2();
     if (a > 0.0) a = 1./std::sqrt(a);
     momentum.setX( a * p.momentum.x() );
     momentum.setY( a * p.momentum.y() );
     momentum.setZ( a * p.momentum.z() );
     mass          = p.mass;
     energy        = std::sqrt(momentum.mag2() + mass*mass);
     kineticEnergy = energy - mass;
     charge        = p.charge;
     timeOfFlight  = p.timeOfFlight;
     side          = p.side;
     flag          = p.flag;
     code          = p.code; 
     particleDef   = p.particleDef; 
     return;
   }

  inline 
   G4double Length()
   {
     return  momentum.mag() ;
   }

  inline
   void Exch( G4GHEKinematicsVector & p1)
   {
     G4GHEKinematicsVector mx = *this;
//     mx.momentum.SetX( momentum.x());
//     mx.momentum.SetY( momentum.y());
//     mx.momentum.SetZ( momentum.z());
//     mx.energy        = energy;
//     mx.kineticEnergy = kineticEnergy;
//     mx.mass          = mass;
//     mx.charge        = charge;
//     mx.timeOfFlight  = timeOfFlight;
//     mx.side          = side;
//     mx.flag          = flag;
//     mx.code          = code; 
//     momentum.setX( p1.momentum.x());
//     momentum.setY( p1.momentum.y());
//     momentum.setZ( p1.momentum.z());
//     energy        = p1.energy;
//     kineticEnergy = p1.kineticEnergy;
//     mass          = p1.mass;
//     charge        = p1.charge;
//     timeOfFlight  = p1.timeOfFlight;
//     side          = p1.side
//     flag          = p1.flag;
//     code          = p1.code;
     *this = p1;
     p1 = mx;
     return; 
   }      

  inline
   void Defs1( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2)
   {
     G4double pt2 = sqr(p1.momentum.x()) + sqr(p1.momentum.y());
     if (pt2 > 0.0)
        {
          G4double ph, px, py, pz;
          G4double cost  = p2.momentum.z()/p2.momentum.mag(); 
          G4double sint  = 0.5 * (  std::sqrt(std::fabs((1.-cost)*(1.+cost))) 
                                  + std::sqrt(pt2)/p2.momentum.mag());
          (p2.momentum.y() < 0.) ? ph = 1.5*CLHEP::pi : ph = CLHEP::halfpi;
          if( p2.momentum.x() != 0.0) 
             ph = std::atan2(p2.momentum.y(),p2.momentum.x());             
          px =   cost*std::cos(ph)*p1.momentum.x() - std::sin(ph)*p1.momentum.y()
               + sint*std::cos(ph)*p1.momentum.z();
          py =   cost*std::sin(ph)*p1.momentum.x() + std::cos(ph)*p1.momentum.y()
               + sint*std::sin(ph)*p1.momentum.z();
          pz = - sint        *p1.momentum.x() 
               + cost        *p1.momentum.z();
          momentum.setX( px ); 
          momentum.setY( py );
          momentum.setZ( pz );     
        }
     else
        {
          momentum = p1.momentum;             
	} 
   }

  inline 
   void Defs( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & p2,
                    G4GHEKinematicsVector & my,       G4GHEKinematicsVector & mz )
   {
     my = p1;
     mz = p2;
     momentum.setX(   my.momentum.y()*mz.momentum.z()
                    - my.momentum.z()*mz.momentum.y());
     momentum.setY(   my.momentum.z()*mz.momentum.x()
                    - my.momentum.x()*mz.momentum.z());
     momentum.setZ(   my.momentum.x()*mz.momentum.y()
                    - my.momentum.y()*mz.momentum.x());
     my.momentum.setX(   mz.momentum.y()*momentum.z()
                       - mz.momentum.z()*momentum.y());
     my.momentum.setY(   mz.momentum.z()*momentum.x()
                       - mz.momentum.x()*momentum.z());
     my.momentum.setZ(   mz.momentum.x()*momentum.y()
                       - mz.momentum.y()*momentum.x());
     G4double pp;
     pp = momentum.mag();
     if (pp > 0.)
        {
          pp = 1./pp; 
          momentum.setX( momentum.x()*pp );
          momentum.setY( momentum.y()*pp );
          momentum.setZ( momentum.z()*pp );
        }
     pp = my.momentum.mag();
     if (pp > 0.)
        {
          pp = 1./pp;
          my.momentum.setX( my.momentum.x()*pp );
          my.momentum.setY( my.momentum.y()*pp );
          my.momentum.setZ( my.momentum.z()*pp );
        }
     pp = mz.momentum.mag();
     if (pp > 0.)
        {
          pp = 1./pp;
          mz.momentum.setX( mz.momentum.x()*pp );
          mz.momentum.setY( mz.momentum.y()*pp );
          mz.momentum.setZ( mz.momentum.z()*pp );
        }
     return; 
   }  
             
  inline
   void Trac( const G4GHEKinematicsVector & p1, const G4GHEKinematicsVector & mx,
              const G4GHEKinematicsVector & my, const G4GHEKinematicsVector & mz)
   {
     double px, py, pz;
     px =   mx.momentum.x()*p1.momentum.x()
          + mx.momentum.y()*p1.momentum.y()
          + mx.momentum.z()*p1.momentum.z();
     py =   my.momentum.x()*p1.momentum.x()
          + my.momentum.y()*p1.momentum.y()
          + my.momentum.z()*p1.momentum.z();
     pz =   mz.momentum.x()*p1.momentum.x()
          + mz.momentum.y()*p1.momentum.y()
          + mz.momentum.z()*p1.momentum.z();
     momentum.setX( px );
     momentum.setY( py );
     momentum.setZ( pz );
     return;
   }                 

  inline
   void Print( G4int LLL)
   {
     G4cout << "G4GHEKinematicsVector: " 
          << LLL << " " << momentum.x() << " " <<  momentum.y() << " " <<  momentum.z() << " "
          << energy << " " << kineticEnergy << " " << mass << " " << charge << " " 
          << timeOfFlight << " " << side << " " << flag << " " << code << particleDef << G4endl;
     return;                         
   }

  G4ParticleMomentum momentum;
  G4double energy;
  G4double kineticEnergy;
  G4double mass;
  G4double charge;
  G4double timeOfFlight;
  G4int side;
  G4bool flag;
  G4int code;
  G4ParticleDefinition * particleDef;
};

#endif                     
