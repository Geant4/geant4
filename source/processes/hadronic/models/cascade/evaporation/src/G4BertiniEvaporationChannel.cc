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
#include "G4BertiniEvaporationChannel.hh"


G4BertiniEvaporationChannel::G4BertiniEvaporationChannel()
{
    verboseLevel = 0;
}


G4BertiniEvaporationChannel::~G4BertiniEvaporationChannel()
{
}


inline void G4BertiniEvaporationChannel::setVerboseLevel( const G4int level )
{
  verboseLevel = level;
}


inline void G4BertiniEvaporationChannel::setNucleusA( const G4int a )
{
  nucleusA = a;
}


inline void G4BertiniEvaporationChannel::setNucleusZ( const G4int z )
{
  nucleusZ = z;
}


inline G4int G4BertiniEvaporationChannel::getParticleA()
{
  return particleA;
}


inline G4int G4BertiniEvaporationChannel::getParticleZ()
{
  return particleZ;
}


inline void G4BertiniEvaporationChannel::setExcitationEnergy(const G4double energy )
{
  excitationEnergy = energy;
}


void G4BertiniEvaporationChannel::setPairingCorrection( const G4int isCorrection )
{
  G4int residualZ = nucleusZ - particleZ;
  G4int residualA = nucleusA - particleA;

  // Note: pairing energy is zero for odd nucleon number.

  if ( isCorrection )  correction = pairingEnergyProtons( residualZ )
	 + pairingEnergyNeutrons( residualA - residualZ );
  else correction = 0;
}


inline G4double G4BertiniEvaporationChannel::getProbability()
{
  if ( verboseLevel >= 4 )
    G4cout << "G4BerEvaporationChannel : probability of particle " << name 
	   << " is " << emissionProbability << endl;

  return emissionProbability;
}


inline G4String G4BertiniEvaporationChannel::getName()
{
  return name;
}


inline void G4BertiniEvaporationChannel::setProbability( const G4double newProb )
{
  emissionProbability = newProb;
  return;
}


G4double G4BertiniEvaporationChannel::getQ()
{
  // The reaction Q is calculated by using the excess masses of the
  // secondary particle and the initial and residual nuclei
  G4int residualZ = nucleusZ  - particleZ;
  G4int residualA = nucleusA  - particleA;
  G4int residualN = residualA - residualZ;
  G4int  nucleusN = nucleusA  - nucleusZ;

  // All the following are by default in MeV
  const G4double e1 = G4NucleiProperties::GetMassExcess( residualA, residualZ );
  const G4double e2 = G4NucleiProperties::GetMassExcess( nucleusA,  nucleusZ  );

  return e1 - e2 + G4NucleiProperties::GetMassExcess( particleA, particleZ ); 
  // In HETC88 : Mass excesses were calculated with the Cameron mass excess formula,
  // see Cameron, Canadian Journal of Physics,
  // vol. 35, 1957, p.1022
  // Also mass excesses of particles were 
  // 8.3675, 7.5851, 13.7279, 15.8381, 15.8195, 3.6092 for
  // n, p, D, T, 3He, 4He
}


G4double G4BertiniEvaporationChannel::getThresh()
{
  G4double threshold = getQ() + getCoulomb();
  return threshold;
}


G4double G4BertiniEvaporationChannel::getCoulomb()
{
  G4int residualZ = nucleusZ - particleZ;
  G4int residualA = nucleusA - particleA;
  const G4double factor = 0.84696; // = e / ( 4 pi epsilon_0 r0 ) * 10^-6, r0=1.7E-15
  // In HETC88 this factor was 0.88235, perhaps due to different r0

  G4double coulomb = factor *  particleZ * qmFactor() * residualZ / 
         ( pow( residualA, 0.33333333 ) + rho ) * MeV;
  
  if ( verboseLevel >= 10 )
    G4cout << " G4BertiniEvaporationChannel::getThresh() " << endl
	   << "        residualZ " << residualZ << endl
	   << "        residualA " << residualA << endl
	   << "         qmFactor " << qmFactor() << endl
	   << "               Q  " << getQ() << endl      
	   << "             rho  " << rho << endl
	   << "           part Z " << particleZ << endl
	   << "     (correction) " << correction << endl;

  return coulomb;
}


inline G4double G4BertiniEvaporationChannel::qmFactor()
{
  //  Coefficient k_p for empirical cross section
  //  formula presented in Dostrovsky, Phys. Rev.,
  //  vol. 116, 1959
   return 0;
}


G4double G4BertiniEvaporationChannel::getLevelDensityParameter()
{
  G4int residualZ = nucleusZ - particleZ;
  G4int residualA = nucleusA - particleA;
  G4double b0 = 8;
  G4double y0 = 1.5;

  G4double temp = ( residualA - 2.0 * residualZ ) / residualA;
  G4double smallA =  residualA * ( 1.0 + y0 * pow( temp, 2 ) ) / b0 / MeV; 

  // In HETC98 b0 = b0(E).

  return smallA;
}


void G4BertiniEvaporationChannel::isotropicCosines( G4double & u, G4double & v, G4double & w )
{
  // Samples isotropic random direction cosines.
  G4double CosTheta = 1.0 - 2.0 * G4UniformRand();
  G4double SinTheta = sqrt( 1.0 - CosTheta * CosTheta );
  G4double Phi = twopi * G4UniformRand();

  u = cos( Phi ) * SinTheta;
  v = cos( Phi ) * CosTheta,
  w = sin( Phi );

  return;
}


G4double G4BertiniEvaporationChannel::pairingEnergyProtons(G4int Z)
{
  //  Pairing energy for protons P(Z), see 
  //  Cameron, Nuclear Level Spacings, p. 1040
  //  Canadian Journal of Physics, vol 36, 1958
  G4double table [130] = {
    0.00000000E+00, 0.54399996E+01, 0.00000000E+00, 0.27599993E+01, 0.00000000E+00,
    0.33399992E+01, 0.00000000E+00, 0.26999998E+01, 0.00000000E+00, 0.18999996E+01,
    0.00000000E+00, 0.21199999E+01, 0.00000000E+00, 0.21299992E+01, 0.00000000E+00,
    0.15399990E+01, 0.00000000E+00, 0.14199991E+01, 0.00000000E+00, 0.15099993E+01,
    0.00000000E+00, 0.17299995E+01, 0.00000000E+00, 0.14399996E+01, 0.00000000E+00,
    0.14499998E+01, 0.00000000E+00, 0.13699999E+01, 0.00000000E+00, 0.10899992E+01,
    0.00000000E+00, 0.13599997E+01, 0.00000000E+00, 0.14199991E+01, 0.00000000E+00,
    0.13299990E+01, 0.00000000E+00, 0.11999998E+01, 0.00000000E+00, 0.99999988E+00,
    0.00000000E+00, 0.11599998E+01, 0.00000000E+00, 0.12799997E+01, 0.00000000E+00,
    0.13799992E+01, 0.00000000E+00, 0.13799992E+01, 0.00000000E+00, 0.13199997E+01,
    0.00000000E+00, 0.10399990E+01, 0.00000000E+00, 0.11099997E+01, 0.00000000E+00,
    0.11299992E+01, 0.00000000E+00, 0.12099991E+01, 0.00000000E+00, 0.14299994E+01,
    0.00000000E+00, 0.11499996E+01, 0.00000000E+00, 0.98999995E+00, 0.00000000E+00,
    0.90999997E+00, 0.00000000E+00, 0.91999996E+00, 0.00000000E+00, 0.99999988E+00,
    0.00000000E+00, 0.11099997E+01, 0.00000000E+00, 0.12299995E+01, 0.00000000E+00,
    0.84999996E+00, 0.00000000E+00, 0.97999996E+00, 0.00000000E+00, 0.71999997E+00,
    0.00000000E+00, 0.79999995E+00, 0.00000000E+00, 0.76999998E+00, 0.00000000E+00,
    0.88999999E+00, 0.00000000E+00, 0.91999996E+00, 0.00000000E+00, 0.79999995E+00,
    0.00000000E+00, 0.80999994E+00, 0.00000000E+00, 0.69000000E+00, 0.00000000E+00,
    0.69999999E+00, 0.00000000E+00, 0.75999999E+00, 0.00000000E+00, 0.72999996E+00,
    0.00000000E+00, 0.79999995E+00, 0.00000000E+00, 0.73999995E+00, 0.00000000E+00,
    0.72999996E+00, 0.00000000E+00, 0.71999997E+00, 0.00000000E+00, 0.71999997E+00,
    0.00000000E+00, 0.71999997E+00, 0.00000000E+00, 0.70999998E+00, 0.00000000E+00,
    0.69000000E+00, 0.00000000E+00, 0.67999995E+00, 0.00000000E+00, 0.65999997E+00,
    0.00000000E+00, 0.60999995E+00, 0.00000000E+00, 0.41999996E+00, 0.00000000E+00,
    0.35999995E+00, 0.00000000E+00, 0.40999997E+00, 0.00000000E+00, 0.48999995E+00
  };
  if ( Z>130 ) G4Exception( " G4BertiniEvaporationChannel: pairing energy for protons called with too large Z " ); 
  return table[ Z-1 ]*MeV;
}


G4double G4BertiniEvaporationChannel::pairingEnergyNeutrons(G4int N)
{
//  Pairing energy for neutrons P(N), see 
//  Cameron, Nuclear Level Spacings, p. 1040
//  Canadian Journal of Physics, vol 36, 1958
  G4double table[200] = {
    0.00000000E+00, 0.59799995E+01, 0.00000000E+00, 0.27699995E+01, 0.00000000E+00,
    0.31599998E+01, 0.00000000E+00, 0.30099993E+01, 0.00000000E+00, 0.16799994E+01,
    0.00000000E+00, 0.17299995E+01, 0.00000000E+00, 0.21699991E+01, 0.00000000E+00,
    0.17399998E+01, 0.00000000E+00, 0.17500000E+01, 0.00000000E+00, 0.17199993E+01,
    0.00000000E+00, 0.16299992E+01, 0.00000000E+00, 0.14099998E+01, 0.00000000E+00,
    0.12899990E+01, 0.00000000E+00, 0.14699993E+01, 0.00000000E+00, 0.13199997E+01,
    0.00000000E+00, 0.14599991E+01, 0.00000000E+00, 0.14399996E+01, 0.00000000E+00,
    0.14599991E+01, 0.00000000E+00, 0.15199995E+01, 0.00000000E+00, 0.15099993E+01,
    0.00000000E+00, 0.14699993E+01, 0.00000000E+00, 0.14499998E+01, 0.00000000E+00,
    0.12799997E+01, 0.00000000E+00, 0.12299995E+01, 0.00000000E+00, 0.12699995E+01,
    0.00000000E+00, 0.61999995E+00, 0.00000000E+00, 0.75999999E+00, 0.00000000E+00,
    0.12299995E+01, 0.00000000E+00, 0.12199993E+01, 0.00000000E+00, 0.13999996E+01,
    0.00000000E+00, 0.13599997E+01, 0.00000000E+00, 0.12999992E+01, 0.00000000E+00,
    0.12899990E+01, 0.00000000E+00, 0.12399998E+01, 0.00000000E+00, 0.12799997E+01,
    0.00000000E+00, 0.12399998E+01, 0.00000000E+00, 0.11999998E+01, 0.00000000E+00,
    0.94000000E+00, 0.00000000E+00, 0.99999988E+00, 0.00000000E+00, 0.10499992E+01,
    0.00000000E+00, 0.53999996E+00, 0.00000000E+00, 0.59999996E+00, 0.00000000E+00,
    0.75000000E+00, 0.00000000E+00, 0.75000000E+00, 0.00000000E+00, 0.84999996E+00,
    0.00000000E+00, 0.96999997E+00, 0.00000000E+00, 0.10199995E+01, 0.00000000E+00,
    0.10499992E+01, 0.00000000E+00, 0.10599995E+01, 0.00000000E+00, 0.10699997E+01,
    0.00000000E+00, 0.10599995E+01, 0.00000000E+00, 0.10499992E+01, 0.00000000E+00,
    0.10199995E+01, 0.00000000E+00, 0.96999997E+00, 0.00000000E+00, 0.90999997E+00,
    0.00000000E+00, 0.82999998E+00, 0.00000000E+00, 0.73999995E+00, 0.00000000E+00,
    0.65999997E+00, 0.00000000E+00, 0.60999995E+00, 0.00000000E+00, 0.60999995E+00,
    0.00000000E+00, 0.89999998E+00, 0.00000000E+00, 0.51999998E+00, 0.00000000E+00,
    0.80999994E+00, 0.00000000E+00, 0.67999995E+00, 0.00000000E+00, 0.71999997E+00,
    0.00000000E+00, 0.76999998E+00, 0.00000000E+00, 0.67999995E+00, 0.00000000E+00,
    0.66999996E+00, 0.00000000E+00, 0.79999995E+00, 0.00000000E+00, 0.67999995E+00,
    0.00000000E+00, 0.63999999E+00, 0.00000000E+00, 0.57999998E+00, 0.00000000E+00,
    0.54999995E+00, 0.00000000E+00, 0.56999993E+00, 0.00000000E+00, 0.56999993E+00,
    0.00000000E+00, 0.54999995E+00, 0.00000000E+00, 0.59999996E+00, 0.00000000E+00,
    0.57999998E+00, 0.00000000E+00, 0.57999998E+00, 0.00000000E+00, 0.60999995E+00,
    0.00000000E+00, 0.63000000E+00, 0.00000000E+00, 0.64999998E+00, 0.00000000E+00,
    0.65999997E+00, 0.00000000E+00, 0.64999998E+00, 0.00000000E+00, 0.64999998E+00,
    0.00000000E+00, 0.63999999E+00, 0.00000000E+00, 0.63999999E+00, 0.00000000E+00,
    0.63000000E+00, 0.00000000E+00, 0.60999995E+00, 0.00000000E+00, 0.58999997E+00,
    0.00000000E+00, 0.54999995E+00, 0.00000000E+00, 0.38999999E+00, 0.00000000E+00,
    0.35999995E+00, 0.00000000E+00, 0.39999998E+00, 0.00000000E+00, 0.39999998E+00,
    0.00000000E+00, 0.39999998E+00, 0.00000000E+00, 0.39999998E+00, 0.00000000E+00,
    0.39999998E+00, 0.00000000E+00, 0.39999998E+00, 0.00000000E+00, 0.39999998E+00
  };
  if ( N > 200 ) G4Exception( " G4BertiniEvaporationChannel: pairing energy for neutrons called with too large Z " ); 
  return table[ N-1 ]*MeV;
}


G4double G4BertiniEvaporationChannel::cameronShellCorrectionP(G4int Z)
{
  // Gives the binding energy correction depending in Z
  // due to shell correction and pairing energies in MeV
  //
  // see Cameron, Canadian Journal of Physics,
  // vol. 35, 1957, p.1022
  G4double table[130] = {
    0.26169998E+02,  0.19250000E+02,  0.24209991E+02,  0.20919998E+02,  0.23149994E+02,
    0.18009995E+02,  0.19549988E+02,  0.16939987E+02,  0.19729996E+02,  0.17069992E+02,
    0.18209991E+02,  0.14990000E+02,  0.16009995E+02,  0.12040000E+02,  0.13270000E+02,
    0.11089998E+02,  0.12169999E+02,  0.10259998E+02,  0.11040000E+02,  0.84099998E+01,
    0.97899990E+01,  0.73599997E+01,  0.81499996E+01,  0.56299992E+01,  0.58799992E+01,
    0.31699991E+01,  0.33199997E+01,  0.81999993E+00,  0.18299999E+01,  0.96999997E+00,
    0.23299999E+01,  0.12699995E+01,  0.29199991E+01,  0.16099997E+01,  0.29099998E+01,
    0.13499994E+01,  0.23999996E+01,  0.88999999E+00,  0.17399998E+01,  0.35999995E+00,
    0.94999999E+00, -0.64999998E+00, -0.39999995E-01, -0.17299995E+01, -0.95999998E+00,
    -0.28699999E+01, -0.20499992E+01, -0.40499992E+01, -0.33999996E+01, -0.57199993E+01,
    -0.37499990E+01, -0.41299992E+01, -0.24199991E+01, -0.28499994E+01, -0.10099993E+01,
    -0.13299990E+01,  0.53999996E+00, -0.20000000E-01,  0.17399998E+01,  0.75000000E+00,
    0.22399998E+01,  0.99999988E+00,  0.19799995E+01,  0.78999996E+00,  0.15399990E+01,
    0.38999999E+00,  0.10799999E+01,  0.00000000E+00,  0.77999997E+00, -0.34999996E+00,
    0.57999998E+00, -0.54999995E+00,  0.58999997E+00, -0.60999995E+00,  0.58999997E+00,
    -0.34999996E+00,  0.31999993E+00, -0.95999998E+00, -0.51999998E+00, -0.20799999E+01,
    -0.24599991E+01, -0.36399994E+01, -0.15499992E+01, -0.95999998E+00,  0.96999997E+00,
    0.88000000E+00,  0.23699999E+01,  0.17500000E+01,  0.27199993E+01,  0.18999996E+01,
    0.25499992E+01,  0.14599991E+01,  0.19299994E+01,  0.85999995E+00,  0.11699991E+01,
    0.79999983E-01,  0.38999999E+00, -0.75999999E+00, -0.38999999E+00, -0.15099993E+01,
    -0.11699991E+01, -0.23599997E+01, -0.19499998E+01, -0.30599995E+01, -0.26199999E+01,
    -0.35499992E+01, -0.29499998E+01, -0.37499990E+01, -0.30699997E+01, -0.37899990E+01,
    -0.30599995E+01, -0.37699995E+01, -0.30499992E+01, -0.37799997E+01, -0.31199999E+01,
    -0.38999996E+01, -0.33499994E+01, -0.42399998E+01, -0.38599997E+01, -0.49199991E+01,
    -0.50599995E+01, -0.67699995E+01, -0.74099998E+01, -0.91799994E+01, -0.10160000E+02,
    -0.11120000E+02, -0.97599993E+01, -0.92299995E+01, -0.79599991E+01, -0.76499996E+01,
  };
  if ( Z > 130 ) G4Exception( " G4BertiniEvaporationChannel: shell correction for protons called with too large Z " ); 
  return table[ Z-1 ]*MeV;
}


G4double G4BertiniEvaporationChannel::cameronShellCorrectionN(G4int N)
{
  // Gives the binding energy correction depending in N
  // due to shell correction and pairing energies in MeV
  //
  // see Cameron, Canadian Journal of Physics,
  // vol. 35, 1957, p.1022
  G4double table[200] = {
    -0.83199997E+01, -0.15899999E+02, -0.11509999E+02, -0.14309999E+02, -0.11570000E+02,
    -0.15899999E+02, -0.13909999E+02, -0.16029999E+02, -0.12129999E+02, -0.13869999E+02,
    -0.12249998E+02, -0.14399999E+02, -0.13069999E+02, -0.15799998E+02, -0.13809999E+02,
    -0.14980000E+02, -0.12629999E+02, -0.13759999E+02, -0.11369999E+02, -0.12379998E+02,
    -0.92299995E+01, -0.96499996E+01, -0.76399994E+01, -0.91699991E+01, -0.80499992E+01,
    -0.97199993E+01, -0.88699989E+01, -0.10759999E+02, -0.86399994E+01, -0.88899994E+01,
    -0.65999994E+01, -0.71299992E+01, -0.47699995E+01, -0.53299990E+01, -0.30599995E+01,
    -0.37899990E+01, -0.17199993E+01, -0.27899990E+01, -0.92999995E+00, -0.21899996E+01,
    -0.51999998E+00, -0.18999996E+01, -0.44999999E+00, -0.21999998E+01, -0.12199993E+01,
    -0.30699997E+01, -0.24199991E+01, -0.43699999E+01, -0.39399996E+01, -0.60799999E+01,
    -0.44899998E+01, -0.45000000E+01, -0.31399994E+01, -0.29299994E+01, -0.10399990E+01,
    -0.13599997E+01,  0.69000000E+00,  0.20999998E+00,  0.21099997E+01,  0.13299990E+01,
    0.32900000E+01,  0.24599991E+01,  0.42999992E+01,  0.33199997E+01,  0.47899990E+01,
    0.36199999E+01,  0.49699993E+01,  0.36399994E+01,  0.46299992E+01,  0.30699997E+01,
    0.40599995E+01,  0.24899998E+01,  0.32999992E+01,  0.14599991E+01,  0.20599995E+01,
    0.50999999E+00,  0.73999995E+00, -0.11799994E+01, -0.12599993E+01, -0.35399990E+01,
    -0.39699993E+01, -0.52599993E+01, -0.41799994E+01, -0.37099991E+01, -0.20999994E+01,
    -0.16999998E+01, -0.79999983E-01, -0.17999995E+00,  0.94000000E+00,  0.26999998E+00,
    0.11299992E+01,  0.79999983E-01,  0.90999997E+00, -0.30999994E+00,  0.48999995E+00,
    -0.77999997E+00,  0.79999983E-01, -0.11499996E+01, -0.22999996E+00, -0.14099998E+01,
    -0.41999996E+00, -0.15499992E+01, -0.54999995E+00, -0.16599998E+01, -0.65999997E+00,
    -0.17299995E+01, -0.75000000E+00, -0.17399998E+01, -0.77999997E+00, -0.16899996E+01,
    -0.77999997E+00, -0.15999994E+01, -0.75000000E+00, -0.14599991E+01, -0.66999996E+00,
    -0.12599993E+01, -0.50999999E+00, -0.10399990E+01, -0.52999997E+00, -0.18399992E+01,
    -0.24199991E+01, -0.45199995E+01, -0.47599993E+01, -0.63299990E+01, -0.67599993E+01,
    -0.78099995E+01, -0.57999992E+01, -0.53699999E+01, -0.36299992E+01, -0.33499994E+01,
    -0.17500000E+01, -0.18799992E+01, -0.60999995E+00, -0.89999998E+00,  0.89999974E-01,
    -0.31999993E+00,  0.54999995E+00, -0.13000000E+00,  0.69999999E+00, -0.59999999E-01,
    0.48999995E+00, -0.19999999E+00,  0.39999998E+00, -0.21999997E+00,  0.35999995E+00,
    -0.89999974E-01,  0.57999998E+00,  0.11999995E+00,  0.75000000E+00,  0.14999998E+00,
    0.69999999E+00,  0.16999996E+00,  0.11099997E+01,  0.88999999E+00,  0.18499994E+01,
    0.16199999E+01,  0.25399990E+01,  0.22899990E+01,  0.31999998E+01,  0.29099998E+01,
    0.38399992E+01,  0.35299997E+01,  0.44799995E+01,  0.41499996E+01,  0.51199999E+01,
    0.47799997E+01,  0.57500000E+01,  0.53899994E+01,  0.63099995E+01,  0.59099998E+01,
    0.68699999E+01,  0.63299990E+01,  0.71299992E+01,  0.66099997E+01,  0.72999992E+01,
    0.63099995E+01,  0.62699995E+01,  0.48299999E+01,  0.44899998E+01,  0.28499994E+01,
    0.23199997E+01,  0.57999998E+00, -0.10999995E+00, -0.97999996E+00,  0.80999994E+00,
    0.17699995E+01,  0.33699999E+01,  0.41299992E+01,  0.55999994E+01,  0.61499996E+01,
    0.72899990E+01,  0.73499994E+01,  0.79499998E+01,  0.76699991E+01,  0.81599998E+01,
    0.78299990E+01,  0.83099995E+01,  0.80099993E+01,  0.85299997E+01,  0.82699995E+01
  };
  if ( N > 130 ) G4Exception( " G4BertiniEvaporationChannel: shell correction for protons called with too large N " ); 
  return table[ N-1 ]*MeV;
}

