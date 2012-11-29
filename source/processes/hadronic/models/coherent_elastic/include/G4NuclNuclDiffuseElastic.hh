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
// $Id$
//
//
// G4 Model: optical elastic scattering with 4-momentum balance 
//
// Class Description
// Final state production model for nucleus-nucleus elastic scattering;
// Coulomb amplitude is not considered as correction 
// (as in G4DiffuseElastic)
// Class Description - End
//
//
// 17.03.09 V. Grichine implementation for Coulomb elastic scattering


#ifndef G4NuclNuclDiffuseElastic_h
#define G4NuclNuclDiffuseElastic_h 1
 
#include <complex>
#include <CLHEP/Units/PhysicalConstants.h>
#include "globals.hh"
#include "G4Integrator.hh"
#include "G4HadronElastic.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

using namespace std;

class G4ParticleDefinition;
class G4PhysicsTable;
class G4PhysicsLogVector;

class G4NuclNuclDiffuseElastic : public G4HadronElastic   // G4HadronicInteraction
{
public:

  G4NuclNuclDiffuseElastic();

  // G4NuclNuclDiffuseElastic(const G4ParticleDefinition* aParticle);





  virtual ~G4NuclNuclDiffuseElastic();

  void Initialise();

  void InitialiseOnFly(G4double Z, G4double A);

  void BuildAngleTable();

 
  // G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, G4Nucleus & targetNucleus);

  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);

  void SetPlabLowLimit(G4double value);

  void SetHEModelLowLimit(G4double value);

  void SetQModelLowLimit(G4double value);

  void SetLowestEnergyLimit(G4double value);

  void SetRecoilKinEnergyLimit(G4double value);

  G4double SampleT(const G4ParticleDefinition* aParticle, 
                         G4double p, G4double A);

  G4double SampleTableT(const G4ParticleDefinition* aParticle, 
                         G4double p, G4double Z, G4double A);

  G4double SampleThetaCMS(const G4ParticleDefinition* aParticle, G4double p, G4double A);

  G4double SampleTableThetaCMS(const G4ParticleDefinition* aParticle, G4double p, 
                                     G4double Z, G4double A);

  G4double GetScatteringAngle(G4int iMomentum, G4int iAngle, G4double position);

  G4double SampleThetaLab(const G4HadProjectile* aParticle, 
                                G4double tmass, G4double A);

  G4double GetDiffuseElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A         );

  G4double GetInvElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A, G4double Z );

  G4double GetDiffuseElasticSumXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A, G4double Z );

  G4double GetInvElasticSumXsc( const G4ParticleDefinition* particle, 
                                 G4double tMand, 
			         G4double momentum, 
				 G4double A, G4double Z );

  G4double IntegralElasticProb( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A            );
  

  G4double GetCoulombElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double Z         );

  G4double GetRutherfordXsc(     G4double theta       );

  G4double GetInvCoulombElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double tMand, 
			         G4double momentum, 
				 G4double A, G4double Z         );

  G4double GetCoulombTotalXsc( const G4ParticleDefinition* particle,  
			         G4double momentum, G4double Z       );

  G4double GetCoulombIntegralXsc( const G4ParticleDefinition* particle,  
			         G4double momentum, G4double Z, 
                                 G4double theta1, G4double theta2         );


  G4double CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    );

  G4double CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 );

  G4double CalculateAm( G4double momentum, G4double n, G4double Z);

  G4double CalculateNuclearRad( G4double A);

  G4double ThetaCMStoThetaLab(const G4DynamicParticle* aParticle, 
                                G4double tmass, G4double thetaCMS);

  G4double ThetaLabToThetaCMS(const G4DynamicParticle* aParticle, 
                                G4double tmass, G4double thetaLab);

  void TestAngleTable(const G4ParticleDefinition* theParticle, G4double partMom,
		      G4double Z, G4double A);



  G4double BesselJzero(G4double z);
  G4double BesselJone(G4double z);
  G4double DampFactor(G4double z);
  G4double BesselOneByArg(G4double z);

  G4double GetDiffElasticProb(G4double theta);
  G4double GetDiffElasticSumProb(G4double theta);
  G4double GetDiffElasticSumProbA(G4double alpha);
  G4double GetIntegrandFunction(G4double theta);

  G4double GetNuclearRadius(){return fNuclearRadius;};


  // Technical math functions for strong Coulomb contribution

  G4complex GammaLogarithm(G4complex xx);
  G4complex GammaLogB2n(G4complex xx);

  G4double  GetErf(G4double x);

  G4double GetCosHaPit2(G4double t){return std::cos(CLHEP::halfpi*t*t);};
  G4double GetSinHaPit2(G4double t){return std::sin(CLHEP::halfpi*t*t);};

  G4double  GetCint(G4double x);
  G4double  GetSint(G4double x);


  G4complex GetErfcComp(G4complex z, G4int nMax);
  G4complex GetErfcSer(G4complex z, G4int nMax);
  G4complex GetErfcInt(G4complex z); // , G4int nMax);

  G4complex GetErfComp(G4complex z, G4int nMax);  // AandS algorithm != Ser, Int
  G4complex GetErfSer(G4complex z, G4int nMax);

  G4double GetExpCos(G4double x);
  G4double GetExpSin(G4double x);
  G4complex GetErfInt(G4complex z); // , G4int nMax);


  

  G4double GetLegendrePol(G4int n, G4double x);

  G4complex TestErfcComp(G4complex z, G4int nMax);
  G4complex TestErfcSer(G4complex z, G4int nMax);
  G4complex TestErfcInt(G4complex z); // , G4int nMax);

  G4complex CoulombAmplitude(G4double theta);
  G4double  CoulombAmplitudeMod2(G4double theta);


  void CalculateCoulombPhaseZero();
  G4double CalculateCoulombPhase(G4int n);
  void CalculateRutherfordAnglePar();

  G4double ProfileNear(G4double theta);
  G4double ProfileFar(G4double theta);
  G4double Profile(G4double theta);

  G4complex PhaseNear(G4double theta);
  G4complex PhaseFar(G4double theta);

  G4complex GammaLess(G4double theta);
  G4complex GammaMore(G4double theta);

  G4complex AmplitudeNear(G4double theta);
  G4complex AmplitudeFar(G4double theta);

  G4complex Amplitude(G4double theta);
  G4double  AmplitudeMod2(G4double theta);

  G4complex AmplitudeSim(G4double theta);
  G4double  AmplitudeSimMod2(G4double theta);

  G4double  GetRatioSim(G4double theta);
  G4double  GetRatioGen(G4double theta);
  
  G4double  GetFresnelDiffuseXsc(G4double theta);
  G4double  GetFresnelIntegrandXsc(G4double alpha);
  

  G4complex AmplitudeGla(G4double theta);
  G4double  AmplitudeGlaMod2(G4double theta);

  G4complex AmplitudeGG(G4double theta);
  G4double  AmplitudeGGMod2(G4double theta);

  void      InitParameters(const G4ParticleDefinition* theParticle,  
			      G4double partMom, G4double Z, G4double A);
 
  void      InitDynParameters(const G4ParticleDefinition* theParticle,  
			      G4double partMom); 

  void      InitParametersGla(const G4DynamicParticle* aParticle,  
			      G4double partMom, G4double Z, G4double A);

  G4double GetHadronNucleonXscNS( G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
				  G4ParticleDefinition* tParticle);

  G4double CalcMandelstamS( const G4double mp , 
                                                       const G4double mt , 
						    const G4double Plab );

  G4double GetProfileLambda(){return fProfileLambda;};

  void SetProfileLambda(G4double pl) {fProfileLambda = pl;};
  void SetProfileDelta(G4double pd) {fProfileDelta = pd;};
  void SetProfileAlpha(G4double pa){fProfileAlpha = pa;};
  void SetCofLambda(G4double pa){fCofLambda = pa;};

  void SetCofAlpha(G4double pa){fCofAlpha = pa;};
  void SetCofAlphaMax(G4double pa){fCofAlphaMax = pa;};
  void SetCofAlphaCoulomb(G4double pa){fCofAlphaCoulomb = pa;};

  void SetCofDelta(G4double pa){fCofDelta = pa;};
  void SetCofPhase(G4double pa){fCofPhase = pa;};
  void SetCofFar(G4double pa){fCofFar = pa;};
  void SetEtaRatio(G4double pa){fEtaRatio = pa;};
  void SetMaxL(G4int l){fMaxL = l;};
  void SetNuclearRadiusCof(G4double r){fNuclearRadiusCof = r;};

  G4double GetCofAlphaMax(){return fCofAlphaMax;};
  G4double GetCofAlphaCoulomb(){return fCofAlphaCoulomb;};

private:


  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theDeuteron;
  G4ParticleDefinition* theAlpha;

  const G4ParticleDefinition* thePionPlus;
  const G4ParticleDefinition* thePionMinus;

  G4double lowEnergyRecoilLimit;  
  G4double lowEnergyLimitHE;  
  G4double lowEnergyLimitQ;  
  G4double lowestEnergyLimit;  
  G4double plabLowLimit;

  G4int fEnergyBin;
  G4int fAngleBin;

  G4PhysicsLogVector*           fEnergyVector;
  G4PhysicsTable*               fAngleTable;
  std::vector<G4PhysicsTable*>  fAngleBank;

  std::vector<G4double> fElementNumberVector;
  std::vector<G4String> fElementNameVector;

  const G4ParticleDefinition* fParticle;

  G4double fWaveVector;
  G4double fAtomicWeight;
  G4double fAtomicNumber;

  G4double fNuclearRadius1;
  G4double fNuclearRadius2;
  G4double fNuclearRadius;
  G4double fNuclearRadiusSquare;
  G4double fNuclearRadiusCof;

  G4double fBeta;
  G4double fZommerfeld;
  G4double fRutherfordRatio;
  G4double fAm;
  G4bool   fAddCoulomb;

  G4double fCoulombPhase0;
  G4double fHalfRutThetaTg;
  G4double fHalfRutThetaTg2;
  G4double fRutherfordTheta;

  G4double fProfileLambda;
  G4double fProfileDelta;
  G4double fProfileAlpha;

  G4double fCofLambda;
  G4double fCofAlpha;
  G4double fCofDelta;
  G4double fCofPhase;
  G4double fCofFar;

  G4double fCofAlphaMax;
  G4double fCofAlphaCoulomb;

  G4int    fMaxL;
  G4double fSumSigma;
  G4double fEtaRatio;

  G4double fReZ;

};


inline void G4NuclNuclDiffuseElastic::SetRecoilKinEnergyLimit(G4double value)
{
  lowEnergyRecoilLimit = value;
}

inline void G4NuclNuclDiffuseElastic::SetPlabLowLimit(G4double value)
{
  plabLowLimit = value;
}

inline void G4NuclNuclDiffuseElastic::SetHEModelLowLimit(G4double value)
{
  lowEnergyLimitHE = value;
}

inline void G4NuclNuclDiffuseElastic::SetQModelLowLimit(G4double value)
{
  lowEnergyLimitQ = value;
}

inline void G4NuclNuclDiffuseElastic::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}


/////////////////////////////////////////////////////////////
//
// Bessel J0 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

inline G4double G4NuclNuclDiffuseElastic::BesselJzero(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = fabs(value);

  if ( value < 8.0 && value > -8.0 )
  {
    value2 = value*value;

    fact1  = 57568490574.0 + value2*(-13362590354.0 
                           + value2*( 651619640.7 
                           + value2*(-11214424.18 
                           + value2*( 77392.33017 
                           + value2*(-184.9052456   ) ) ) ) );

    fact2  = 57568490411.0 + value2*( 1029532985.0 
                           + value2*( 9494680.718
                           + value2*(59272.64853
                           + value2*(267.8532712 
                           + value2*1.0               ) ) ) );

    bessel = fact1/fact2;
  } 
  else 
  {
    arg    = 8.0/modvalue;

    value2 = arg*arg;

    shift  = modvalue-0.785398164;

    fact1  = 1.0 + value2*(-0.1098628627e-2 
                 + value2*(0.2734510407e-4
                 + value2*(-0.2073370639e-5 
                 + value2*0.2093887211e-6    ) ) );

    fact2  = -0.1562499995e-1 + value2*(0.1430488765e-3
                              + value2*(-0.6911147651e-5 
                              + value2*(0.7621095161e-6
                              - value2*0.934945152e-7    ) ) );

    bessel = sqrt(0.636619772/modvalue)*(cos(shift)*fact1 - arg*sin(shift)*fact2 );
  }
  return bessel;
}

/////////////////////////////////////////////////////////////
//
// Bessel J1 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

inline G4double G4NuclNuclDiffuseElastic::BesselJone(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = fabs(value);

  if ( modvalue < 8.0 ) 
  {
    value2 = value*value;

    fact1  = value*(72362614232.0 + value2*(-7895059235.0 
                                  + value2*( 242396853.1
                                  + value2*(-2972611.439 
                                  + value2*( 15704.48260 
                                  + value2*(-30.16036606  ) ) ) ) ) );

    fact2  = 144725228442.0 + value2*(2300535178.0 
                            + value2*(18583304.74
                            + value2*(99447.43394 
                            + value2*(376.9991397 
                            + value2*1.0             ) ) ) );
    bessel = fact1/fact2;
  } 
  else 
  {
    arg    = 8.0/modvalue;

    value2 = arg*arg;

    shift  = modvalue - 2.356194491;

    fact1  = 1.0 + value2*( 0.183105e-2 
                 + value2*(-0.3516396496e-4
                 + value2*(0.2457520174e-5 
                 + value2*(-0.240337019e-6          ) ) ) );

    fact2 = 0.04687499995 + value2*(-0.2002690873e-3
                          + value2*( 0.8449199096e-5
                          + value2*(-0.88228987e-6
                          + value2*0.105787412e-6       ) ) );

    bessel = sqrt( 0.636619772/modvalue)*(cos(shift)*fact1 - arg*sin(shift)*fact2);

    if (value < 0.0) bessel = -bessel;
  }
  return bessel;
}

////////////////////////////////////////////////////////////////////
//
// damp factor in diffraction x/sh(x), x was already *pi

inline G4double G4NuclNuclDiffuseElastic::DampFactor(G4double x)
{
  G4double df;
  G4double f2 = 2., f3 = 6., f4 = 24.; // first factorials

  // x *= pi;

  if( std::fabs(x) < 0.01 )
  { 
    df = 1./(1. + x/f2 + x*x/f3 + x*x*x/f4);
  }
  else
  {
    df = x/std::sinh(x); 
  }
  return df;
}


////////////////////////////////////////////////////////////////////
//
// return J1(x)/x with special case for small x

inline G4double G4NuclNuclDiffuseElastic::BesselOneByArg(G4double x)
{
  G4double x2, result;
  
  if( std::fabs(x) < 0.01 )
  { 
   x     *= 0.5;
   x2     = x*x;
   result = 2. - x2 + x2*x2/6.;
  }
  else
  {
    result = BesselJone(x)/x; 
  }
  return result;
}

////////////////////////////////////////////////////////////////////
//
// return particle beta

inline  G4double G4NuclNuclDiffuseElastic::CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    )
{
  G4double mass = particle->GetPDGMass();
  G4double a    = momentum/mass;
  fBeta         = a/std::sqrt(1+a*a);

  return fBeta; 
}

////////////////////////////////////////////////////////////////////
//
// return Zommerfeld parameter for Coulomb scattering

inline  G4double G4NuclNuclDiffuseElastic::CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 )
{
  fZommerfeld = CLHEP::fine_structure_const*Z1*Z2/beta;

  return fZommerfeld; 
}

////////////////////////////////////////////////////////////////////
//
// return Wentzel correction for Coulomb scattering

inline  G4double G4NuclNuclDiffuseElastic::CalculateAm( G4double momentum, G4double n, G4double Z)
{
  G4double k   = momentum/CLHEP::hbarc;
  G4double ch  = 1.13 + 3.76*n*n;
  G4double zn  = 1.77*k*std::pow(Z,-1./3.)*CLHEP::Bohr_radius;
  G4double zn2 = zn*zn;
  fAm          = ch/zn2;

  return fAm;
}

////////////////////////////////////////////////////////////////////
//
// calculate nuclear radius for different atomic weights using different approximations

inline  G4double G4NuclNuclDiffuseElastic::CalculateNuclearRad( G4double A)
{
  G4double r0 = 1.*CLHEP::fermi, radius;
  // r0 *= 1.12;
  // r0 *= 1.44;
  r0 *= fNuclearRadiusCof;

  /*
  if( A < 50. )
  {
    if( A > 10. ) r0  = 1.16*( 1 - std::pow(A, -2./3.) )*CLHEP::fermi;   // 1.08*fermi;
    else          r0  = 1.1*CLHEP::fermi;

    radius = r0*std::pow(A, 1./3.);
  }
  else
  {
    r0 = 1.7*CLHEP::fermi;   // 1.7*fermi;

    radius = r0*std::pow(A, 0.27); // 0.27);
  }
  */
  radius = r0*std::pow(A, 1./3.);

  return radius;
}

////////////////////////////////////////////////////////////////////
//
// return Coulomb scattering differential xsc with Wentzel correction. Test function  

inline  G4double G4NuclNuclDiffuseElastic::GetCoulombElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double Z         )
{
  G4double sinHalfTheta  = std::sin(0.5*theta);
  G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta;
  G4double beta          = CalculateParticleBeta( particle, momentum);
  G4double z             = particle->GetPDGCharge();
  G4double n             = CalculateZommerfeld( beta, z, Z );
  G4double am            = CalculateAm( momentum, n, Z);
  G4double k             = momentum/CLHEP::hbarc;
  G4double ch            = 0.5*n/k;
  G4double ch2           = ch*ch;
  G4double xsc           = ch2/(sinHalfTheta2+am)/(sinHalfTheta2+am);

  return xsc;
}

////////////////////////////////////////////////////////////////////
//
// return Rutherford scattering differential xsc with Wentzel correction. For Sampling.  

inline  G4double G4NuclNuclDiffuseElastic::GetRutherfordXsc(   G4double theta  )
{
  G4double sinHalfTheta  = std::sin(0.5*theta);
  G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta;

  G4double ch2           = fRutherfordRatio*fRutherfordRatio;

  G4double xsc           = ch2/(sinHalfTheta2+fAm)/(sinHalfTheta2+fAm);

  return xsc;
}


////////////////////////////////////////////////////////////////////
//
// return Coulomb scattering total xsc with Wentzel correction  

inline  G4double G4NuclNuclDiffuseElastic::GetCoulombTotalXsc( const G4ParticleDefinition* particle,  
			                                     G4double momentum, G4double Z  )
{
  G4double beta          = CalculateParticleBeta( particle, momentum);
  G4cout<<"beta = "<<beta<<G4endl;
  G4double z             = particle->GetPDGCharge();
  G4double n             = CalculateZommerfeld( beta, z, Z );
  G4cout<<"fZomerfeld = "<<n<<G4endl;
  G4double am            = CalculateAm( momentum, n, Z);
  G4cout<<"cof Am = "<<am<<G4endl;
  G4double k             = momentum/CLHEP::hbarc;
  G4cout<<"k = "<<k*CLHEP::fermi<<" 1/fermi"<<G4endl;
  G4cout<<"k*Bohr_radius = "<<k*CLHEP::Bohr_radius<<G4endl;
  G4double ch            = n/k;
  G4double ch2           = ch*ch;
  G4double xsc           = ch2*CLHEP::pi/(am +am*am);

  return xsc;
}

////////////////////////////////////////////////////////////////////
//
// return Coulomb scattering xsc with Wentzel correction  integrated between
// theta1 and < theta2

inline  G4double G4NuclNuclDiffuseElastic::GetCoulombIntegralXsc( const G4ParticleDefinition* particle,  
			         G4double momentum, G4double Z, 
                                 G4double theta1, G4double theta2 )
{
  G4double c1 = std::cos(theta1);
  G4cout<<"c1 = "<<c1<<G4endl;
  G4double c2 = std::cos(theta2);
  G4cout<<"c2 = "<<c2<<G4endl;
  G4double beta          = CalculateParticleBeta( particle, momentum);
  // G4cout<<"beta = "<<beta<<G4endl;
  G4double z             = particle->GetPDGCharge();
  G4double n             = CalculateZommerfeld( beta, z, Z );
  // G4cout<<"fZomerfeld = "<<n<<G4endl;
  G4double am            = CalculateAm( momentum, n, Z);
  // G4cout<<"cof Am = "<<am<<G4endl;
  G4double k             = momentum/CLHEP::hbarc;
  // G4cout<<"k = "<<k*CLHEP::fermi<<" 1/fermi"<<G4endl;
  // G4cout<<"k*Bohr_radius = "<<k*CLHEP::Bohr_radius<<G4endl;
  G4double ch            = n/k;
  G4double ch2           = ch*ch;
  am *= 2.;
  G4double xsc           = ch2*CLHEP::twopi*(c1-c2);
           xsc          /= (1 - c1 + am)*(1 - c2 + am);

  return xsc;
}

///////////////////////////////////////////////////////////////////
//
// For the calculation of arg Gamma(z) one needs complex extension 
// of ln(Gamma(z))

inline G4complex G4NuclNuclDiffuseElastic::GammaLogarithm(G4complex zz)
{
  static G4double cof[6] = { 76.18009172947146,     -86.50532032941677,
                             24.01409824083091,      -1.231739572450155,
                              0.1208650973866179e-2, -0.5395239384953e-5  } ;
  register G4int j;
  G4complex z = zz - 1.0;
  G4complex tmp = z + 5.5;
  tmp -= (z + 0.5) * std::log(tmp);
  G4complex ser = G4complex(1.000000000190015,0.);

  for ( j = 0; j <= 5; j++ )
  {
    z += 1.0;
    ser += cof[j]/z;
  }
  return -tmp + std::log(2.5066282746310005*ser);
}

///////////////////////////////////////////////////////////////////
//
// For the calculation of arg Gamma(z) one needs complex extension 
// of ln(Gamma(z)) here is approximate algorithm

inline G4complex G4NuclNuclDiffuseElastic::GammaLogB2n(G4complex z)
{
  G4complex z1 = 12.*z;
  G4complex z2 = z*z;
  G4complex z3 = z2*z;
  G4complex z5 = z2*z3;
  G4complex z7 = z2*z5;

  z3 *= 360.;
  z5 *= 1260.;
  z7 *= 1680.;

  G4complex result  = (z-0.5)*std::log(z)-z+0.5*std::log(CLHEP::twopi);
            result += 1./z1 - 1./z3 +1./z5 -1./z7;
  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4double  G4NuclNuclDiffuseElastic::GetErf(G4double x)
{
  G4double t, z, tmp, result;

  z   = std::fabs(x);
  t   = 1.0/(1.0+0.5*z);

  tmp = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));

  if( x >= 0.) result = 1. - tmp;
  else         result = 1. + tmp; 
    
  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfcComp(G4complex z, G4int nMax)
{
  G4complex erfcz = 1. - GetErfComp( z, nMax);
  return erfcz;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfcSer(G4complex z, G4int nMax)
{
  G4complex erfcz = 1. - GetErfSer( z, nMax);
  return erfcz;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfcInt(G4complex z) // , G4int nMax)
{
  G4complex erfcz = 1. - GetErfInt( z); // , nMax);
  return erfcz;
}

inline  G4double G4NuclNuclDiffuseElastic::GetLegendrePol(G4int n, G4double theta)
{
  G4double legPol, epsilon = 1.e-6;
  G4double x = std::cos(theta);

  if     ( n  < 0 ) legPol = 0.;
  else if( n == 0 ) legPol = 1.;
  else if( n == 1 ) legPol = x;
  else if( n == 2 ) legPol = (3.*x*x-1.)/2.;
  else if( n == 3 ) legPol = (5.*x*x*x-3.*x)/2.;
  else if( n == 4 ) legPol = (35.*x*x*x*x-30.*x*x+3.)/8.;
  else if( n == 5 ) legPol = (63.*x*x*x*x*x-70.*x*x*x+15.*x)/8.;
  else if( n == 6 ) legPol = (231.*x*x*x*x*x*x-315.*x*x*x*x+105.*x*x-5.)/16.;
  else           
  {
    // legPol = ( (2*n-1)*x*GetLegendrePol(n-1,x) - (n-1)*GetLegendrePol(n-2,x) )/n;

    legPol = std::sqrt( 2./(n*CLHEP::pi*std::sin(theta+epsilon)) )*std::sin( (n+0.5)*theta+0.25*CLHEP::pi );
  }
  return legPol; 
}



/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::TestErfcComp(G4complex z, G4int nMax)
{
  G4complex miz = G4complex( z.imag(), -z.real() ); 
  G4complex erfcz = 1. - GetErfComp( miz, nMax);
  G4complex w = std::exp(-z*z)*erfcz;
  return w;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::TestErfcSer(G4complex z, G4int nMax)
{
  G4complex miz = G4complex( z.imag(), -z.real() ); 
  G4complex erfcz = 1. - GetErfSer( miz, nMax);
  G4complex w = std::exp(-z*z)*erfcz;
  return w;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::TestErfcInt(G4complex z) // , G4int nMax)
{
  G4complex miz = G4complex( z.imag(), -z.real() ); 
  G4complex erfcz = 1. - GetErfInt( miz); // , nMax);
  G4complex w = std::exp(-z*z)*erfcz;
  return w;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfComp(G4complex z, G4int nMax)
{
  G4int n;
  G4double n2, cofn, shny, chny, fn, gn;

  G4double x = z.real();
  G4double y = z.imag();

  G4double outRe = 0., outIm = 0.;

  G4double twox  = 2.*x;
  G4double twoxy = twox*y;
  G4double twox2 = twox*twox;

  G4double cof1 = std::exp(-x*x)/CLHEP::pi;

  G4double cos2xy = std::cos(twoxy);
  G4double sin2xy = std::sin(twoxy);

  G4double twoxcos2xy = twox*cos2xy;
  G4double twoxsin2xy = twox*sin2xy;

  for( n = 1; n <= nMax; n++)
  {
    n2   = n*n;

    cofn = std::exp(-0.5*n2)/(n2+twox2);  // /(n2+0.5*twox2);

    chny = std::cosh(n*y);
    shny = std::sinh(n*y);

    fn  = twox - twoxcos2xy*chny + n*sin2xy*shny;
    gn  =        twoxsin2xy*chny + n*cos2xy*shny;

    fn *= cofn;
    gn *= cofn;

    outRe += fn;
    outIm += gn;
  }
  outRe *= 2*cof1;
  outIm *= 2*cof1;

  if(std::abs(x) < 0.0001)
  {
    outRe += GetErf(x);
    outIm += cof1*y;
  }
  else
  {
    outRe += GetErf(x) + cof1*(1-cos2xy)/twox;
    outIm += cof1*sin2xy/twox;
  }
  return G4complex(outRe, outIm);
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfSer(G4complex z, G4int nMax)
{
  G4int n;
  G4double a =1., b = 1., tmp;
  G4complex sum = z, d = z;

  for( n = 1; n <= nMax; n++)
  {
    a *= 2.;
    b *= 2.*n +1.;
    d *= z*z;

    tmp = a/b;

    sum += tmp*d;
  }
  sum *= 2.*std::exp(-z*z)/std::sqrt(CLHEP::pi);

  return sum;
}

/////////////////////////////////////////////////////////////////////

inline  G4double  G4NuclNuclDiffuseElastic::GetExpCos(G4double x)
{
  G4double result;

  result = std::exp(x*x-fReZ*fReZ);
  result *= std::cos(2.*x*fReZ);
  return result;
}

/////////////////////////////////////////////////////////////////////

inline  G4double  G4NuclNuclDiffuseElastic::GetExpSin(G4double x)
{
  G4double result;

  result = std::exp(x*x-fReZ*fReZ);
  result *= std::sin(2.*x*fReZ);
  return result;
}



/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GetErfInt(G4complex z) // , G4int nMax)
{
  G4double outRe, outIm;

  G4double x = z.real();
  G4double y = z.imag();
  fReZ       = x;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  outRe = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetExpSin, 0., y );
  outIm = integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetExpCos, 0., y );

  outRe *= 2./sqrt(CLHEP::pi);
  outIm *= 2./sqrt(CLHEP::pi);

  outRe += GetErf(x);

  return G4complex(outRe, outIm);
}

/////////////////////////////////////////////////////////////////
//
//

inline G4double G4NuclNuclDiffuseElastic::GetCint(G4double x)
{
  G4double out;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  out= integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetCosHaPit2, 0., x );

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4double G4NuclNuclDiffuseElastic::GetSint(G4double x)
{
  G4double out;

  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  out= integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetSinHaPit2, 0., x );

  return out;
}


/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::CoulombAmplitude(G4double theta)
{
  G4complex ca;
  
  G4double sinHalfTheta  = std::sin(0.5*theta);
  G4double sinHalfTheta2 = sinHalfTheta*sinHalfTheta; 
  sinHalfTheta2         += fAm;

  G4double order         = 2.*fCoulombPhase0 - fZommerfeld*std::log(sinHalfTheta2);
  G4complex z            = G4complex(0., order);
  ca                     = std::exp(z);

  ca                    *= -fZommerfeld/(2.*fWaveVector*sinHalfTheta2);

  return ca; 
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4double G4NuclNuclDiffuseElastic::CoulombAmplitudeMod2(G4double theta)
{
  G4complex ca = CoulombAmplitude(theta);
  G4double out = ca.real()*ca.real() + ca.imag()*ca.imag();

  return  out;
}

/////////////////////////////////////////////////////////////////
//
//


inline  void G4NuclNuclDiffuseElastic::CalculateCoulombPhaseZero()
{
  G4complex z        = G4complex(1,fZommerfeld); 
  // G4complex gammalog = GammaLogarithm(z);
  G4complex gammalog = GammaLogB2n(z);
  fCoulombPhase0     = gammalog.imag();
}

/////////////////////////////////////////////////////////////////
//
//


inline  G4double G4NuclNuclDiffuseElastic::CalculateCoulombPhase(G4int n)
{
  G4complex z        = G4complex(1. + n, fZommerfeld); 
  // G4complex gammalog = GammaLogarithm(z);
  G4complex gammalog = GammaLogB2n(z);
  return gammalog.imag();
}


/////////////////////////////////////////////////////////////////
//
//


inline  void G4NuclNuclDiffuseElastic::CalculateRutherfordAnglePar()
{
  fHalfRutThetaTg   = fZommerfeld/fProfileLambda;  // (fWaveVector*fNuclearRadius);
  fRutherfordTheta  = 2.*std::atan(fHalfRutThetaTg);
  fHalfRutThetaTg2  = fHalfRutThetaTg*fHalfRutThetaTg;
  // G4cout<<"fRutherfordTheta = "<<fRutherfordTheta/degree<<" degree"<<G4endl;

}

/////////////////////////////////////////////////////////////////
//
//

inline   G4double G4NuclNuclDiffuseElastic::ProfileNear(G4double theta)
{
  G4double dTheta = fRutherfordTheta - theta;
  G4double result = 0., argument = 0.;

  if(std::abs(dTheta) < 0.001) result = fProfileAlpha*fProfileDelta;
  else
  {
    argument = fProfileDelta*dTheta;
    result   = CLHEP::pi*argument*std::exp(fProfileAlpha*argument);
    result  /= std::sinh(CLHEP::pi*argument);
    result  -= 1.;
    result  /= dTheta;
  }
  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline   G4double G4NuclNuclDiffuseElastic::ProfileFar(G4double theta)
{
  G4double dTheta   = fRutherfordTheta + theta;
  G4double argument = fProfileDelta*dTheta;

  G4double result   = CLHEP::pi*argument*std::exp(fProfileAlpha*argument);
  result           /= std::sinh(CLHEP::pi*argument);
  result           /= dTheta;

  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline   G4double G4NuclNuclDiffuseElastic::Profile(G4double theta)
{
  G4double dTheta = fRutherfordTheta - theta;
  G4double result = 0., argument = 0.;

  if(std::abs(dTheta) < 0.001) result = 1.;
  else
  {
    argument = fProfileDelta*dTheta;
    result   = CLHEP::pi*argument;
    result  /= std::sinh(CLHEP::pi*argument);
  }
  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline   G4complex G4NuclNuclDiffuseElastic::PhaseNear(G4double theta)
{
  G4double twosigma = 2.*fCoulombPhase0; 
  twosigma -= fZommerfeld*std::log(fHalfRutThetaTg2/(1.+fHalfRutThetaTg2));
  twosigma += fRutherfordTheta*fZommerfeld/fHalfRutThetaTg - CLHEP::halfpi;
  twosigma -= fProfileLambda*theta - 0.25*CLHEP::pi;

  twosigma *= fCofPhase;

  G4complex z = G4complex(0., twosigma);

  return std::exp(z);
}

/////////////////////////////////////////////////////////////////
//
//

inline   G4complex G4NuclNuclDiffuseElastic::PhaseFar(G4double theta)
{
  G4double twosigma = 2.*fCoulombPhase0; 
  twosigma -= fZommerfeld*std::log(fHalfRutThetaTg2/(1.+fHalfRutThetaTg2));
  twosigma += fRutherfordTheta*fZommerfeld/fHalfRutThetaTg - CLHEP::halfpi;
  twosigma += fProfileLambda*theta - 0.25*CLHEP::pi;

  twosigma *= fCofPhase;

  G4complex z = G4complex(0., twosigma);

  return std::exp(z);
}

/////////////////////////////////////////////////////////////////
//
//


inline G4complex G4NuclNuclDiffuseElastic::GammaLess(G4double theta)
{
  G4double sinThetaR      = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double cosHalfThetaR2 = 1./(1. + fHalfRutThetaTg2);

  G4double u              = std::sqrt(0.5*fProfileLambda/sinThetaR);
  G4double kappa          = u/std::sqrt(CLHEP::pi);
  G4double dTheta         = theta - fRutherfordTheta;
  u                      *= dTheta;
  G4double u2             = u*u;
  G4double u2m2p3         = u2*2./3.;

  G4complex im            = G4complex(0.,1.);
  G4complex order         = G4complex(u,u);
  order                  /= std::sqrt(2.);

  G4complex gamma         = CLHEP::pi*kappa*GetErfcInt(-order)*std::exp(im*(u*u+0.25*CLHEP::pi));
  G4complex a0            = 0.5*(1. + 4.*(1.+im*u2)*cosHalfThetaR2/3.)/sinThetaR;
  G4complex a1            = 0.5*(1. + 2.*(1.+im*u2m2p3)*cosHalfThetaR2)/sinThetaR;
  G4complex out           = gamma*(1. - a1*dTheta) - a0;

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::GammaMore(G4double theta)
{
  G4double sinThetaR      = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double cosHalfThetaR2 = 1./(1. + fHalfRutThetaTg2);

  G4double u              = std::sqrt(0.5*fProfileLambda/sinThetaR);
  G4double kappa          = u/std::sqrt(CLHEP::pi);
  G4double dTheta         = theta - fRutherfordTheta;
  u                      *= dTheta;
  G4double u2             = u*u;
  G4double u2m2p3         = u2*2./3.;

  G4complex im            = G4complex(0.,1.);
  G4complex order         = G4complex(u,u);
  order                  /= std::sqrt(2.);
  G4complex gamma         = CLHEP::pi*kappa*GetErfcInt(order)*std::exp(im*(u*u+0.25*CLHEP::pi));
  G4complex a0            = 0.5*(1. + 4.*(1.+im*u2)*cosHalfThetaR2/3.)/sinThetaR;
  G4complex a1            = 0.5*(1. + 2.*(1.+im*u2m2p3)*cosHalfThetaR2)/sinThetaR;
  G4complex out           = -gamma*(1. - a1*dTheta) - a0;

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::AmplitudeNear(G4double theta)
{
  G4double kappa = std::sqrt(0.5*fProfileLambda/std::sin(theta)/CLHEP::pi);
  G4complex out = G4complex(kappa/fWaveVector,0.);

  out *= PhaseNear(theta);

  if( theta <= fRutherfordTheta )
  {
    out *= GammaLess(theta) + ProfileNear(theta);
    // out *= GammaMore(theta) + ProfileNear(theta);
    out += CoulombAmplitude(theta);
  }
  else
  {
    out *= GammaMore(theta) + ProfileNear(theta);
    // out *= GammaLess(theta) + ProfileNear(theta);
  }
  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::AmplitudeFar(G4double theta)
{
  G4double kappa = std::sqrt(0.5*fProfileLambda/std::sin(theta)/CLHEP::pi);
  G4complex out = G4complex(kappa/fWaveVector,0.);
  out *= ProfileFar(theta);
  out *= PhaseFar(theta);
  return out;
}


/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::Amplitude(G4double theta)
{
 
  G4complex out = AmplitudeNear(theta) + fCofFar*AmplitudeFar(theta);
  // G4complex out = AmplitudeNear(theta);
  // G4complex out = AmplitudeFar(theta);
  return    out;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4double  G4NuclNuclDiffuseElastic::AmplitudeMod2(G4double theta)
{
  G4complex out = Amplitude(theta);
  G4double mod2 = out.real()*out.real() + out.imag()*out.imag();
  return   mod2;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4complex G4NuclNuclDiffuseElastic::AmplitudeSim(G4double theta)
{
  G4double sinThetaR  = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double dTheta     = 0.5*(theta - fRutherfordTheta);
  G4double sindTheta  = std::sin(dTheta);
  G4double persqrt2   = std::sqrt(0.5);

  G4complex order     = G4complex(persqrt2,persqrt2);
  order              *= std::sqrt(0.5*fProfileLambda/sinThetaR)*2.*sindTheta;
  // order              *= std::sqrt(0.5*fProfileLambda/sinThetaR)*2.*dTheta;

  G4complex out;

  if ( theta <= fRutherfordTheta )
  {
    out = 1. - 0.5*GetErfcInt(-order)*ProfileNear(theta);
  }
  else
  {
    out = 0.5*GetErfcInt(order)*ProfileNear(theta);
  }

  out *= CoulombAmplitude(theta);

  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline G4double G4NuclNuclDiffuseElastic::GetRatioSim(G4double theta)
{
  G4double sinThetaR  = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double dTheta     = 0.5*(theta - fRutherfordTheta);
  G4double sindTheta  = std::sin(dTheta);

  G4double order      = std::sqrt(fProfileLambda/sinThetaR/CLHEP::pi)*2.*sindTheta;
  // G4cout<<"order = "<<order<<G4endl;  
  G4double cosFresnel = 0.5 - GetCint(order);  
  G4double sinFresnel = 0.5 - GetSint(order);  

  G4double out = 0.5*( cosFresnel*cosFresnel + sinFresnel*sinFresnel );

  return out;
}

/////////////////////////////////////////////////////////////////
//
// The ratio el/ruth for Fresnel smooth nucleus profile

inline G4double G4NuclNuclDiffuseElastic::GetRatioGen(G4double theta)
{
  G4double sinThetaR  = 2.*fHalfRutThetaTg/(1. + fHalfRutThetaTg2);
  G4double dTheta     = 0.5*(theta - fRutherfordTheta);
  G4double sindTheta  = std::sin(dTheta);

  G4double prof       = Profile(theta);
  G4double prof2      = prof*prof;
  // G4double profmod    = std::abs(prof);
  G4double order      = std::sqrt(fProfileLambda/sinThetaR/CLHEP::pi)*2.*sindTheta;

  order = std::abs(order); // since sin changes sign!
  // G4cout<<"order = "<<order<<G4endl; 
 
  G4double cosFresnel = GetCint(order);  
  G4double sinFresnel = GetSint(order);  

  G4double out;

  if ( theta <= fRutherfordTheta )
  {
    out  = 1. + 0.5*( (0.5-cosFresnel)*(0.5-cosFresnel)+(0.5-sinFresnel)*(0.5-sinFresnel) )*prof2; 
    out += ( cosFresnel + sinFresnel - 1. )*prof;
  }
  else
  {
    out = 0.5*( (0.5-cosFresnel)*(0.5-cosFresnel)+(0.5-sinFresnel)*(0.5-sinFresnel) )*prof2;
  }

  return out;
}

/////////////////////////////////////////////////////////////////
//
// The xsc for Fresnel smooth nucleus profile

inline G4double G4NuclNuclDiffuseElastic::GetFresnelDiffuseXsc(G4double theta)
{
  G4double ratio   = GetRatioGen(theta);
  G4double ruthXsc = GetRutherfordXsc(theta);
  G4double xsc     = ratio*ruthXsc;
  return xsc;
}

/////////////////////////////////////////////////////////////////
//
// The xsc for Fresnel smooth nucleus profile for integration

inline G4double G4NuclNuclDiffuseElastic::GetFresnelIntegrandXsc(G4double alpha)
{
  G4double theta = std::sqrt(alpha);
  G4double xsc     = GetFresnelDiffuseXsc(theta);
  return xsc;
}



/////////////////////////////////////////////////////////////////
//
//

inline  G4double  G4NuclNuclDiffuseElastic::AmplitudeSimMod2(G4double theta)
{
  G4complex out = AmplitudeSim(theta);
  G4double mod2 = out.real()*out.real() + out.imag()*out.imag();
  return   mod2;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::AmplitudeGla(G4double theta)
{
  G4int n;
  G4double T12b, b, b2; // cosTheta = std::cos(theta);
  G4complex out = G4complex(0.,0.), shiftC, shiftN; 
  G4complex im  = G4complex(0.,1.);

  for( n = 0; n < fMaxL; n++)
  {
    shiftC = std::exp( im*2.*CalculateCoulombPhase(n) );
    // b = ( fZommerfeld + std::sqrt( fZommerfeld*fZommerfeld + n*(n+1) ) )/fWaveVector;
    b = ( std::sqrt( G4double(n*(n+1)) ) )/fWaveVector;
    b2 = b*b;
    T12b = fSumSigma*std::exp(-b2/fNuclearRadiusSquare)/CLHEP::pi/fNuclearRadiusSquare;         
    shiftN = std::exp( -0.5*(1.-im*fEtaRatio)*T12b ) - 1.;
    out +=  (2.*n+1.)*shiftC*shiftN*GetLegendrePol(n, theta);   
  }
  out /= 2.*im*fWaveVector;
  out += CoulombAmplitude(theta);
  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4double  G4NuclNuclDiffuseElastic::AmplitudeGlaMod2(G4double theta)
{
  G4complex out = AmplitudeGla(theta);
  G4double mod2 = out.real()*out.real() + out.imag()*out.imag();
  return   mod2;
}


/////////////////////////////////////////////////////////////////
//
//

inline  G4complex G4NuclNuclDiffuseElastic::AmplitudeGG(G4double theta)
{
  G4int n;
  G4double T12b, a, aTemp, b2, sinThetaH = std::sin(0.5*theta);
  G4double sinThetaH2 = sinThetaH*sinThetaH;
  G4complex out = G4complex(0.,0.); 
  G4complex im  = G4complex(0.,1.);

  a  = -fSumSigma/CLHEP::twopi/fNuclearRadiusSquare;
  b2 = fWaveVector*fWaveVector*fNuclearRadiusSquare*sinThetaH2;

  aTemp = a;

  for( n = 1; n < fMaxL; n++)
  {    
    T12b   = aTemp*std::exp(-b2/n)/n;         
    aTemp *= a;  
    out   += T12b; 
    G4cout<<"out = "<<out<<G4endl;  
  }
  out *= -4.*im*fWaveVector/CLHEP::pi;
  out += CoulombAmplitude(theta);
  return out;
}

/////////////////////////////////////////////////////////////////
//
//

inline  G4double  G4NuclNuclDiffuseElastic::AmplitudeGGMod2(G4double theta)
{
  G4complex out = AmplitudeGG(theta);
  G4double mod2 = out.real()*out.real() + out.imag()*out.imag();
  return   mod2;
}


///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

inline void G4NuclNuclDiffuseElastic::InitParameters(const G4ParticleDefinition* theParticle,  
                                          G4double partMom, G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // atomic number
  fAtomicWeight  = A;     // number of nucleons

  fNuclearRadius2 = CalculateNuclearRad(fAtomicWeight);
  G4double A1     = G4double( theParticle->GetBaryonNumber() );   
  fNuclearRadius1 = CalculateNuclearRad(A1);
  // fNuclearRadius = std::sqrt(fNuclearRadius1*fNuclearRadius1+fNuclearRadius2*fNuclearRadius2);
  fNuclearRadius = fNuclearRadius1 + fNuclearRadius2;

  G4double a = 0.;
  G4double z = theParticle->GetPDGCharge();
  G4double m1 = theParticle->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double lambda = fCofLambda*fWaveVector*fNuclearRadius;
  G4cout<<"kR = "<<lambda<<G4endl;

  if( z )
  {
    a           = partMom/m1; // beta*gamma for m1
    fBeta       = a/std::sqrt(1+a*a);
    fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
    fRutherfordRatio = fZommerfeld/fWaveVector; 
    fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }
  G4cout<<"fZommerfeld = "<<fZommerfeld<<G4endl;
  fProfileLambda = lambda; // *std::sqrt(1.-2*fZommerfeld/lambda);
  G4cout<<"fProfileLambda = "<<fProfileLambda<<G4endl;
  fProfileDelta  = fCofDelta*fProfileLambda;
  fProfileAlpha  = fCofAlpha*fProfileLambda;

  CalculateCoulombPhaseZero();
  CalculateRutherfordAnglePar();

  return;
}
///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

inline void G4NuclNuclDiffuseElastic::InitDynParameters(const G4ParticleDefinition* theParticle,  
                                          G4double partMom) 
{
  G4double a = 0.;
  G4double z = theParticle->GetPDGCharge();
  G4double m1 = theParticle->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double lambda = fCofLambda*fWaveVector*fNuclearRadius;

  if( z )
  {
    a           = partMom/m1; // beta*gamma for m1
    fBeta       = a/std::sqrt(1+a*a);
    fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
    fRutherfordRatio = fZommerfeld/fWaveVector; 
    fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }
  fProfileLambda = lambda; // *std::sqrt(1.-2*fZommerfeld/lambda);
  fProfileDelta  = fCofDelta*fProfileLambda;
  fProfileAlpha  = fCofAlpha*fProfileLambda;

  CalculateCoulombPhaseZero();
  CalculateRutherfordAnglePar();

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
// Test for given particle and element table of momentum, angle probability.
// For the partMom in CMS. 

inline void G4NuclNuclDiffuseElastic::InitParametersGla(const G4DynamicParticle* aParticle,  
                                          G4double partMom, G4double Z, G4double A) 
{
  fAtomicNumber  = Z;     // target atomic number
  fAtomicWeight  = A;     // target number of nucleons

  fNuclearRadius2 = CalculateNuclearRad(fAtomicWeight); // target nucleus radius
  G4double A1     = G4double( aParticle->GetDefinition()->GetBaryonNumber() );   
  fNuclearRadius1 = CalculateNuclearRad(A1); // projectile nucleus radius
  fNuclearRadiusSquare = fNuclearRadius1*fNuclearRadius1+fNuclearRadius2*fNuclearRadius2;
 

  G4double a  = 0., kR12;
  G4double z  = aParticle->GetDefinition()->GetPDGCharge();
  G4double m1 = aParticle->GetDefinition()->GetPDGMass();

  fWaveVector = partMom/CLHEP::hbarc;

  G4double pN = A1 - z;
  if( pN < 0. ) pN = 0.;

  G4double tN = A - Z;
  if( tN < 0. ) tN = 0.;

  G4double pTkin = aParticle->GetKineticEnergy();  
  pTkin /= A1;


  fSumSigma = (Z*z+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
              (z*tN+pN*Z)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  G4cout<<"fSumSigma = "<<fSumSigma/CLHEP::millibarn<<" mb"<<G4endl;
  G4cout<<"pi*R2 = "<<CLHEP::pi*fNuclearRadiusSquare/CLHEP::millibarn<<" mb"<<G4endl;
  kR12 = fWaveVector*std::sqrt(fNuclearRadiusSquare);
  G4cout<<"k*sqrt(R2) = "<<kR12<<" "<<G4endl;
  fMaxL = (G4int(kR12)+1)*4;
  G4cout<<"fMaxL = "<<fMaxL<<" "<<G4endl;

  if( z )
  {
      a           = partMom/m1; // beta*gamma for m1
      fBeta       = a/std::sqrt(1+a*a);
      fZommerfeld = CalculateZommerfeld( fBeta, z, fAtomicNumber);
      fAm         = CalculateAm( partMom, fZommerfeld, fAtomicNumber);
  }

  CalculateCoulombPhaseZero();
 

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns nucleon-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
// projectile nucleon is pParticle with pTkin shooting target nucleon tParticle

inline G4double 
G4NuclNuclDiffuseElastic::GetHadronNucleonXscNS( G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
                                                 G4ParticleDefinition* tParticle)
{
  G4double xsection(0), /*Delta,*/ A0, B0;
  G4double hpXsc(0);
  G4double hnXsc(0);


  G4double targ_mass     = tParticle->GetPDGMass();
  G4double proj_mass     = pParticle->GetPDGMass(); 

  G4double proj_energy   = proj_mass + pTkin; 
  G4double proj_momentum = std::sqrt(pTkin*(pTkin+2*proj_mass));

  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= CLHEP::GeV*CLHEP::GeV;  // in GeV for parametrisation
  proj_momentum /= CLHEP::GeV;
  proj_energy   /= CLHEP::GeV;
  proj_mass     /= CLHEP::GeV;
  G4double logS = std::log(sMand);

  // General PDG fit constants


  // fEtaRatio=Re[f(0)]/Im[f(0)]

  if( proj_momentum >= 1.2 )
  {
    fEtaRatio  = 0.13*(logS - 5.8579332)*std::pow(sMand,-0.18);
  }       
  else if( proj_momentum >= 0.6 )
  { 
    fEtaRatio = -75.5*(std::pow(proj_momentum,0.25)-0.95)/
	  (std::pow(3*proj_momentum,2.2)+1);     
  }
  else 
  {
    fEtaRatio = 15.5*proj_momentum/(27*proj_momentum*proj_momentum*proj_momentum+2);
  }
  G4cout<<"fEtaRatio = "<<fEtaRatio<<G4endl;

  // xsc
  
  if( proj_momentum >= 10. ) // high energy: pp = nn = np
    // if( proj_momentum >= 2.)
  {
    //Delta = 1.;

    //if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

    if( proj_momentum >= 10.)
    {
        B0 = 7.5;
        A0 = 100. - B0*std::log(3.0e7);

        xsection = A0 + B0*std::log(proj_energy) - 11
                  + 103*std::pow(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
    }
  }
  else // low energy pp = nn != np
  {
      if(pParticle == tParticle) // pp or nn      // nn to be pp
      {
        if( proj_momentum < 0.73 )
        {
          hnXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hnXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
        }
        else  // if( proj_momentum < 10.  )
        {
          hnXsc = 39.0 + 
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
        }
        xsection = hnXsc;
      }
      else  // pn to be np
      {
        if( proj_momentum < 0.8 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
        }
        else    // if( proj_momentum < 10.  )
        {
          hpXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
        }
        xsection = hpXsc;
      }
  }
  xsection *= CLHEP::millibarn; // parametrised in mb
  G4cout<<"xsection = "<<xsection/CLHEP::millibarn<<" mb"<<G4endl;
  return xsection;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

inline G4double G4NuclNuclDiffuseElastic::CalcMandelstamS( const G4double mp , 
                                                       const G4double mt , 
                                                       const G4double Plab )
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}



#endif
