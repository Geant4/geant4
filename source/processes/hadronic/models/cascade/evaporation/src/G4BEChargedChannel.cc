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

#include "G4BEChargedChannel.hh"


G4BEChargedChannel::G4BEChargedChannel()
{
    verboseLevel = 0;
}


G4BEChargedChannel::~G4BEChargedChannel()
{
}


void G4BEChargedChannel::calculateProbability()
{
  G4int residualZ = nucleusZ - particleZ;
  G4int residualA = nucleusA - particleA;

//    Check if nucleus is too small, if this evaporation channel
//    leads to an impossible residual nucleus or if there is no enough
//    energy.
  if ( nucleusA < 2.0 * particleA || 
       nucleusZ < 2.0 * particleZ ||
       residualA <= residualZ || 
       excitationEnergy - getThresh() - correction < 0 )
    {
      if ( verboseLevel >= 6 )
	G4cout << "G4BEChargedChannel : calculateProbability for " << getName() << " = 0 " << endl;
      emissionProbability = 0;
      return;
    }

  // In HETC88 s-s0 was used in exp( s ), in which s0 was either 50 or
  // max(s_i), where i goes over all channels.

  G4double levelParam = getLevelDensityParameter();
  G4double s        = 2 * sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
  G4double constant = A / 2 * ( 2 * spin + 1 ) * ( 1 + coulombFactor() );
  G4double eye1     = ( pow( s, 2 ) - 3 * s + 3 ) / ( 4 * pow( levelParam, 2 ) ) * exp( s );

  emissionProbability = constant * pow( residualA, 0.6666666 ) * eye1;

  if ( verboseLevel >= 6 )
    G4cout << "G4BEChargedChannel : calculateProbability for " << getName() << endl
	   << "                    res A = " << residualA << endl 
	   << "                    res Z = " << residualZ << endl 
	   << "                 c*(c_i+1) = "<< constant << endl
	   << "                  qmfactor = "<< qmFactor() << endl
	   << "             coulombfactor = "<< coulombFactor() << endl
	   << "                        E = " << excitationEnergy << endl
	   << "               correction = " << correction << endl
	   << "                     eye1 = " << eye1 << endl
	   << "               levelParam = " << levelParam << endl
	   << "                   thresh = " << getThresh() << endl
	   << "                        s = " << s << endl
	   << "              probability = " << emissionProbability << endl;

  return;
}


G4double G4BEChargedChannel::sampleKineticEnergy()
{
//    G4double randExp1;
//    G4double randExp2;
//    G4double s;
//    G4double levelParam;
//    G4double kineticEnergyAv;
//    G4double kineticEnergy;
  
//    randExp1 = RandExponential::shoot( 1 );
//    randExp2 = RandExponential::shoot( 1 );
//    levelParam = getLevelDensityParameter();
//    s = 2 * sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
//    kineticEnergyAv = 2 * ( pow( s, 3 ) - 6.0 * pow( s, 2 ) + 15.0 * s - 15.0 )  / 
//        ( ( 2.0 * pow( s, 2 ) - 6.0 * s + 6.0 ) * levelParam );
  
//    kineticEnergy = 0.5 * ( randExp1 + randExp2 ) * kineticEnergyAv + getThresh() - getQ();

//    if ( verboseLevel >= 10 )
//      G4cout << "  G4BEChargedChannel : sampleKineticEnergy() " << endl
//  	   << "         kinetic e = " << kineticEnergy << endl
//  	   << "           average = " << kineticEnergyAv << endl
//  	   << "                 s = " << s << endl
//  	   << "        levelParam = " << levelParam << endl
//  	   << "          randExp1 = " << randExp1 << endl
//  	   << "          randExp2 = " << randExp2 << endl;

  G4double levelParam;
  levelParam = getLevelDensityParameter();
  
  const G4double xMax  = excitationEnergy - getThresh() - correction; // maximum number
  const G4double xProb = ( - 1 + sqrt ( 1 + 4 * levelParam * xMax ) ) / ( 2 * levelParam ); // most probable value
  const G4double m = xProb * exp ( 2 * sqrt ( levelParam * ( xMax - xProb ) ) ); // maximum value of P(x)

  // Sample x according to density function P(x) with rejection method
  G4double r1;
  G4double r2;
  G4int koe=0;
  do
    {
      r1 = G4UniformRand() * xMax;
      r2 = G4UniformRand() * m;
      koe++;
    }
  while (  r1 * exp ( 2 * sqrt ( levelParam * ( xMax - r1 ) ) )  < r2 );

  G4cout << "Q ch " << koe << endl;
  G4double kineticEnergy = r1 + getCoulomb(); // add coulomb potential;

  if ( verboseLevel >= 10 )
    G4cout << " G4BENeutronChannel : sampleKineticEnergy() " << endl
	   << "       kinetic n e = " << kineticEnergy << endl
	   << "        levelParam = " << levelParam << endl
	   << "             thresh= " << getThresh() << endl;

  return kineticEnergy;
}


G4double G4BEChargedChannel::coulombFactorForProton()
{
  //  Coefficient c_p:s for empirical cross section formula are
  //  defined with the proton constant.  See Dostrovsky, Phys. Rev.,
  //  vol. 116, 1959.
  G4double t[7] = { 0.08 , 0 , -0.06 , -0.1 , -0.1 , -0.1 , -0.1 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];
  
  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  G4double ret_val =  x * t[n - 2] + ( 1.0 - x ) * t[n-1];

  return ret_val;
}


G4double G4BEChargedChannel::qmFactorForProton()
{
  //  Coefficient k_p for empirical cross section formula are defined
  //  with the proton constant.  See Dostrovsky, Phys. Rev., vol. 116,
  //  1959
  G4double t[7] = {  0.36, 0.51, 0.60, 0.66, 0.68, 0.69, 0.69 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];

  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  return x * t[n - 2] + ( 1.0 - x ) * t[n-1];
}


G4double G4BEChargedChannel::qmFactorForAlpha()
{
//  Coefficient k_alpha for empirical cross section formula presented
//  in Dostrovsky, Phys. Rev., vol. 116, 1959

  G4double t[7] = {  0.77, 0.81, 0.85, 0.89, 0.93, 0.97,  1.00 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];

  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  return x * t[n - 2] + ( 1.0 - x ) * t[n-1];
}
