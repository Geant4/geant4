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
// $Id: G4NuclNuclDiffuseElastic.hh 106722 2017-10-20 09:48:19Z gcosmo $
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

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"


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

  G4double SampleThetaCMS( const G4ParticleDefinition* aParticle, G4double p, G4double A);

  G4double SampleCoulombMuCMS( const G4ParticleDefinition* aParticle, G4double p);

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

  G4double CalcMandelstamS( const G4double mp , const G4double mt , 
			    const G4double Plab );

  G4double GetProfileLambda(){return fProfileLambda;};

  inline void SetProfileLambda(G4double pl) {fProfileLambda = pl;};
  inline void SetProfileDelta(G4double pd) {fProfileDelta = pd;};
  inline void SetProfileAlpha(G4double pa){fProfileAlpha = pa;};
  inline void SetCofLambda(G4double pa){fCofLambda = pa;};

  inline void SetCofAlpha(G4double pa){fCofAlpha = pa;};
  inline void SetCofAlphaMax(G4double pa){fCofAlphaMax = pa;};
  inline void SetCofAlphaCoulomb(G4double pa){fCofAlphaCoulomb = pa;};

  inline void SetCofDelta(G4double pa){fCofDelta = pa;};
  inline void SetCofPhase(G4double pa){fCofPhase = pa;};
  inline void SetCofFar(G4double pa){fCofFar = pa;};
  inline void SetEtaRatio(G4double pa){fEtaRatio = pa;};
  inline void SetMaxL(G4int l){fMaxL = l;};
  inline void SetNuclearRadiusCof(G4double r){fNuclearRadiusCof = r;};

  inline G4double GetCofAlphaMax(){return fCofAlphaMax;};
  inline G4double GetCofAlphaCoulomb(){return fCofAlphaCoulomb;};

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
  G4double fCoulombMuC;

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

inline  G4double 
G4NuclNuclDiffuseElastic::CalculateParticleBeta(const G4ParticleDefinition* particle, 
						G4double momentum)
{
  G4double mass = particle->GetPDGMass();
  G4double a    = momentum/mass;
  fBeta         = a/std::sqrt(1+a*a);

  return fBeta; 
}

////////////////////////////////////////////////////////////////////
//
// return Zommerfeld parameter for Coulomb scattering

inline  G4double 
G4NuclNuclDiffuseElastic::CalculateZommerfeld(G4double beta, G4double Z1, G4double Z2 )
{
  fZommerfeld = CLHEP::fine_structure_const*Z1*Z2/beta;

  return fZommerfeld; 
}

////////////////////////////////////////////////////////////////////
//
// return Wentzel correction for Coulomb scattering

inline  G4double 
G4NuclNuclDiffuseElastic::CalculateAm( G4double momentum, G4double n, G4double Z)
{
  G4double k   = momentum/CLHEP::hbarc;
  G4double ch  = 1.13 + 3.76*n*n;
  G4double zn  = 1.77*k*(1.0/G4Pow::GetInstance()->A13(Z))*CLHEP::Bohr_radius;
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
    if( A > 10. ) r0  = 1.16*( 1 - (1.0/G4Pow::GetInstance()->A23(A)) )*CLHEP::fermi;   // 1.08*fermi;
    else          r0  = 1.1*CLHEP::fermi;

    radius = r0*G4Pow::GetInstance()->A13(A);
  }
  else
  {
    r0 = 1.7*CLHEP::fermi;   // 1.7*fermi;

    radius = r0*G4Pow::GetInstance()->powA(A, 0.27); // 0.27);
  }
  */
  radius = r0*G4Pow::GetInstance()->A13(A);

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
  G4double xsc           = ch2/((sinHalfTheta2+am)*(sinHalfTheta2+am));

  return xsc;
}

////////////////////////////////////////////////////////////////////
//
// return Rutherford scattering differential xsc with Wentzel correction. For Sampling.  

inline  G4double G4NuclNuclDiffuseElastic::GetRutherfordXsc(G4double theta)
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

inline  G4double 
G4NuclNuclDiffuseElastic::GetCoulombIntegralXsc(const G4ParticleDefinition* particle,  
						G4double momentum, G4double Z, 
						G4double theta1, G4double theta2 )
{
  G4double c1 = std::cos(theta1);
  //G4cout<<"c1 = "<<c1<<G4endl;
  G4double c2 = std::cos(theta2);
  // G4cout<<"c2 = "<<c2<<G4endl;
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
  G4double xsc = ch2*CLHEP::twopi*(c1-c2)/((1 - c1 + am)*(1 - c2 + am));

  return xsc;
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

  G4complex result  = (z-0.5)*std::log(z)-z+0.5*G4Log(CLHEP::twopi);
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

  tmp = t*std::exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
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

  result = G4Exp(x*x-fReZ*fReZ);
  result *= std::cos(2.*x*fReZ);
  return result;
}

/////////////////////////////////////////////////////////////////////

inline  G4double  G4NuclNuclDiffuseElastic::GetExpSin(G4double x)
{
  G4double result;

  result = G4Exp(x*x-fReZ*fReZ);
  result *= std::sin(2.*x*fReZ);
  return result;
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
  G4Integrator<G4NuclNuclDiffuseElastic,G4double(G4NuclNuclDiffuseElastic::*)(G4double)> integral;

  G4double out = 
    integral.Legendre96(this,&G4NuclNuclDiffuseElastic::GetSinHaPit2, 0., x );

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

  G4double order         = 2.*fCoulombPhase0 - fZommerfeld*G4Log(sinHalfTheta2);
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
    result   = CLHEP::pi*argument*G4Exp(fProfileAlpha*argument);
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

  G4double result   = CLHEP::pi*argument*G4Exp(fProfileAlpha*argument);
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
    result  /= std::sinh(result);
  }
  return result;
}

/////////////////////////////////////////////////////////////////
//
//

inline   G4complex G4NuclNuclDiffuseElastic::PhaseNear(G4double theta)
{
  G4double twosigma = 2.*fCoulombPhase0; 
  twosigma -= fZommerfeld*G4Log(fHalfRutThetaTg2/(1.+fHalfRutThetaTg2));
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
  twosigma -= fZommerfeld*G4Log(fHalfRutThetaTg2/(1.+fHalfRutThetaTg2));
  twosigma += fRutherfordTheta*fZommerfeld/fHalfRutThetaTg - CLHEP::halfpi;
  twosigma += fProfileLambda*theta - 0.25*CLHEP::pi;

  twosigma *= fCofPhase;

  G4complex z = G4complex(0., twosigma);

  return std::exp(z);
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

inline  G4double  G4NuclNuclDiffuseElastic::AmplitudeGlaMod2(G4double theta)
{
  G4complex out = AmplitudeGla(theta);
  G4double mod2 = out.real()*out.real() + out.imag()*out.imag();
  return   mod2;
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
