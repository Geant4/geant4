// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ParticleTable.hh"    
#include "G4Nucleus.hh"  
#include "G4BENeutronChannel.hh"


G4BENeutronChannel::G4BENeutronChannel()
{
  name         = "neutron";
  particleA    = 1;
  particleZ    = 0;
  verboseLevel = 0;
  rho          = 0;
}


G4BENeutronChannel::~G4BENeutronChannel()
{
}


void G4BENeutronChannel::calculateProbability()
{
  const G4int residualZ = nucleusZ - particleZ;
  const G4int residualA = nucleusA - particleA;

  if ( nucleusA  <  2.0 * particleA ||  
       nucleusZ  <  2.0 * particleZ ||  
       residualA <= residualZ       ||  
       excitationEnergy - getThresh() - correction < 0 )
    {
      if ( verboseLevel >= 6 )
	G4cout << "G4BENeutronChannel : calculateProbability = 0 " << endl;
      emissionProbability = 0;
      return;
    }

  // In HETC88 s-s0 was used in exp( s ), in which s0 was either 50 or
  // max(s_i), where i goes over all channels.

  const G4double levelParam = getLevelDensityParameter();
  
  const G4double s    = 2 * sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
  const G4double temp = ( pow( s, 2 ) - 3 * s + 3 ) / ( 4 * pow( levelParam, 2 ) ) 
    + beta() * ( s - 1 ) /  ( 2 * levelParam );
  const G4double eye0 = exp( s ) * ( s - 1 ) /  ( 2 * levelParam );
  const G4double eye1 = ( pow( s, 2 ) - 3*s +3 ) * exp( s ) / ( 4 * pow( levelParam, 2 ) ) ;
  
  emissionProbability = pow( residualA, 0.666666 ) * alpha() * ( eye1 + beta() * eye0 );
  
  if ( verboseLevel >= 6 )
    G4cout << "G4BENeutronChannel : calculateProbability " << endl
	   << "                    res A = " << residualA << endl 
	   << "                    res Z = " << residualZ << endl 
	   << "                    alpha = " << alpha() << endl 
	   << "                     beta = " << beta() << endl
	   << "                        E = " << excitationEnergy << endl
	   << "               correction = " << correction << endl
	   << "                     eye1 = " << eye1 << endl
	   << "                     eye0 = " <<  eye0 << endl
	   << "               levelParam = " << levelParam << endl
	   << "                   thresh = " << getThresh() << endl
	   << "                        s = " << s << endl
	   << "              probability = " << emissionProbability << endl;

  return;
}


G4double  G4BENeutronChannel::sampleKineticEnergy()
{
  // Samples the kinetic energy of the particle in CMS
  //
  // Algorithm used in HETC98
  //
//    G4double e1;
//    G4double e2;
//    G4double s;
//    G4double levelParam;
//    G4double eye0;
//    G4double eye1;
//    G4double kineticEnergyAv;
//    G4double kineticEnergy;
  
//    e1 = RandExponential::shoot( 1 );
//    e2 = RandExponential::shoot( 1 );

//    levelParam  = getLevelDensityParameter();
//    s = 2 * sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
//    eye0 = 0.5 * ( s - 1 ) * exp( s ) / levelParam;
//    eye1 = ( pow( s, 2 ) - 3*s + 3 ) * exp( s ) / ( 4 * pow( levelParam, 2 ) );
//    kineticEnergyAv = 2 * ( pow( s, 3 ) - 6.0 * pow( s, 2 ) + 15.0 * s - 15.0 )  / 
//        ( ( 2.0 * pow( s, 2 ) - 6.0 * s + 6.0 ) * levelParam );
//    kineticEnergyAv = ( kineticEnergyAv + beta() ) / ( 1.0 + beta() * eye0
//  	   / eye1 );
  
//    kineticEnergy = 0.5 * ( e1 + e2 ) * kineticEnergyAv + getThresh() - getQ();

  ////////////////
  // A random number is sampled from the density function
  // P(x) = x * exp ( 2 sqrt ( a ( xMax - x ) ) )  [not normalized],
  // x belongs to [ 0, xMax ]
  // with the 'Hit or Miss' -method
  // Kinetic energy is this energy scaled properly
  
  G4double levelParam;
  levelParam = getLevelDensityParameter();
  
  const G4double xMax  = excitationEnergy - getThresh() - correction + beta(); // maximum number
  const G4double xProb = ( - 1 + sqrt ( 1 + 4 * levelParam * xMax ) ) / ( 2 * levelParam ); // most probable value
  const G4double m = xProb * exp ( 2 * sqrt ( levelParam * ( xMax - xProb ) ) ); // maximum value of P(x)

  // Sample x according to density function P(x) with rejection method
  G4double r1;
  G4double r2;
  G4int koe=0;
  do
    {
      r1 = beta() + G4UniformRand() * ( xMax - beta() );
      r2 = G4UniformRand() * m;
      koe++;
    }
  while (  r1 * exp ( 2 * sqrt ( levelParam * ( xMax - r1 ) ) )  < r2 );

  G4cout <<  koe << endl;
  G4double kineticEnergy = r1 - beta();

  if ( verboseLevel >= 10 )
    G4cout << " G4BENeutronChannel : sampleKineticEnergy() " << endl
	   << "       kinetic n e = " << kineticEnergy << endl
	   << "        levelParam = " << levelParam << endl
	   << "             thresh= " << getThresh() << endl
	   << "               beta= " << beta() << endl;

  return kineticEnergy;
}


G4DynamicParticle *  G4BENeutronChannel::emit()
{
  G4double u;
  G4double v;
  G4double w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;

  pParticle -> SetDefinition( G4Neutron::Neutron() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u , v , w );  

  return pParticle;
}


G4double G4BENeutronChannel::alpha()
{
  const G4double residualA = nucleusA - particleA;
  return 0.76 + 1.93 * pow( residualA, -0.33333 );
}


G4double G4BENeutronChannel::beta()
{
  G4double residualA = nucleusA - particleA;
  return ( 1.66 * pow ( residualA, -0.66666 ) - 0.05 )/alpha()*MeV;
}

