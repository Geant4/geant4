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
// Author: V. Grichine (Vladimir,Grichine@cern.ch)
//
//
// G4 Model: diffuse optical elastic scattering with 4-momentum balance 
//
// Class Description
// Final state production model for hadron nuclear elastic scattering; 
// Class Description - End
//
//
// 24.05.07 V. Grichine, first implementation for hadron (no Coulomb) elastic scattering
// 04.09.07 V. Grichine, implementation for Coulomb elastic scattering
// 12.06.11 V. Grichine, new interface to G4hadronElastic
// 24.11.17 W. Pokorski, code cleanup and performance improvements


#ifndef G4DiffuseElasticV2_h
#define G4DiffuseElasticV2_h 1

#include <CLHEP/Units/PhysicalConstants.h>
#include "globals.hh"
#include "G4HadronElastic.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

#include "G4Pow.hh"

#include <vector>

class G4ParticleDefinition;
class G4PhysicsTable;
class G4PhysicsLogVector;

class G4DiffuseElasticV2 : public G4HadronElastic // G4HadronicInteraction
{
public:

  G4DiffuseElasticV2();

  virtual ~G4DiffuseElasticV2();

  virtual G4bool IsApplicable(const G4HadProjectile &/*aTrack*/, 
			      G4Nucleus & /*targetNucleus*/);

  void Initialise();

  void InitialiseOnFly(G4double Z, G4double A);

  void BuildAngleTable();

  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);

  G4double NeutronTuniform(G4int Z);

  void SetPlabLowLimit(G4double value);

  void SetHEModelLowLimit(G4double value);

  void SetQModelLowLimit(G4double value);

  void SetLowestEnergyLimit(G4double value);

  void SetRecoilKinEnergyLimit(G4double value);

  G4double SampleTableT(const G4ParticleDefinition* aParticle, 
                         G4double p, G4double Z, G4double A);

  G4double SampleThetaCMS(const G4ParticleDefinition* aParticle, G4double p, G4double A);

  G4double SampleTableThetaCMS(const G4ParticleDefinition* aParticle, G4double p, 
                                     G4double Z, G4double A);

  G4double GetScatteringAngle(G4int iMomentum, unsigned long iAngle, G4double position);

  G4double SampleThetaLab(const G4HadProjectile* aParticle, 
                                G4double tmass, G4double A);

  G4double CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 );

  G4double CalculateAm( G4double momentum, G4double n, G4double Z);

  G4double CalculateNuclearRad( G4double A);

  G4double ThetaCMStoThetaLab(const G4DynamicParticle* aParticle, 
                                G4double tmass, G4double thetaCMS);

  G4double ThetaLabToThetaCMS(const G4DynamicParticle* aParticle, 
                                G4double tmass, G4double thetaLab);


  G4double BesselJzero(G4double z);
  G4double BesselJone(G4double z);
  G4double DampFactor(G4double z);
  G4double BesselOneByArg(G4double z);

  G4double GetDiffElasticSumProbA(G4double alpha);
  G4double GetIntegrandFunction(G4double theta);


  G4double GetNuclearRadius(){return fNuclearRadius;};

private:


  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;

  G4double lowEnergyRecoilLimit;  
  G4double lowEnergyLimitHE;  
  G4double lowEnergyLimitQ;  
  G4double lowestEnergyLimit;  
  G4double plabLowLimit;

  G4int fEnergyBin;
  unsigned long fAngleBin;

  G4PhysicsLogVector* fEnergyVector;
  
  std::vector<std::vector<std::vector<double>*>*>  fEnergyAngleVectorBank;
  std::vector<std::vector<std::vector<double>*>*>  fEnergySumVectorBank;

  std::vector<std::vector<double>*>*  fEnergyAngleVector;
  std::vector<std::vector<double>*>*  fEnergySumVector;

  
  std::vector<G4double> fElementNumberVector;
  std::vector<G4String> fElementNameVector;

  const G4ParticleDefinition* fParticle;
  G4double fWaveVector;
  G4double fAtomicWeight;
  G4double fAtomicNumber;
  G4double fNuclearRadius;
  G4double fBeta;
  G4double fZommerfeld;
  G4double fAm;
  G4bool fAddCoulomb;

};

inline G4bool G4DiffuseElasticV2::IsApplicable(const G4HadProjectile & projectile, 
			      G4Nucleus & nucleus)
{
  if( ( projectile.GetDefinition() == G4Proton::Proton() ||
        projectile.GetDefinition() == G4Neutron::Neutron() ||
        projectile.GetDefinition() == G4PionPlus::PionPlus() ||
        projectile.GetDefinition() == G4PionMinus::PionMinus() ||
        projectile.GetDefinition() == G4KaonPlus::KaonPlus() ||
        projectile.GetDefinition() == G4KaonMinus::KaonMinus() ) &&

        nucleus.GetZ_asInt() >= 2 ) return true;
  else                              return false;
}

inline void G4DiffuseElasticV2::SetRecoilKinEnergyLimit(G4double value)
{
  lowEnergyRecoilLimit = value;
}

inline void G4DiffuseElasticV2::SetPlabLowLimit(G4double value)
{
  plabLowLimit = value;
}

inline void G4DiffuseElasticV2::SetHEModelLowLimit(G4double value)
{
  lowEnergyLimitHE = value;
}

inline void G4DiffuseElasticV2::SetQModelLowLimit(G4double value)
{
  lowEnergyLimitQ = value;
}

inline void G4DiffuseElasticV2::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}


/////////////////////////////////////////////////////////////
//
// Bessel J0 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

inline G4double G4DiffuseElasticV2::BesselJzero(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = std::fabs(value);

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

    bessel = std::sqrt(0.636619772/modvalue)*(std::cos(shift)*fact1 - arg*std::sin(shift)*fact2 );
  }
  return bessel;
}

/////////////////////////////////////////////////////////////
//
// Bessel J1 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

inline G4double G4DiffuseElasticV2::BesselJone(G4double value)
{
  G4double modvalue, value2, fact1, fact2, arg, shift, bessel;

  modvalue = std::fabs(value);

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

    bessel = std::sqrt( 0.636619772/modvalue)*(std::cos(shift)*fact1 - arg*std::sin(shift)*fact2);

    if (value < 0.0) bessel = -bessel;
  }
  return bessel;
}

////////////////////////////////////////////////////////////////////
//
// damp factor in diffraction x/sh(x), x was already *pi

inline G4double G4DiffuseElasticV2::DampFactor(G4double x)
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

inline G4double G4DiffuseElasticV2::BesselOneByArg(G4double x)
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
// return Zommerfeld parameter for Coulomb scattering

inline  G4double G4DiffuseElasticV2::CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 )
{
  fZommerfeld = CLHEP::fine_structure_const*Z1*Z2/beta;

  return fZommerfeld; 
}

////////////////////////////////////////////////////////////////////
//
// return Wentzel correction for Coulomb scattering

inline  G4double G4DiffuseElasticV2::CalculateAm( G4double momentum, G4double n, G4double Z)
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

inline  G4double G4DiffuseElasticV2::CalculateNuclearRad( G4double A)
{
  G4double R, r0, a11, a12, a13, a2, a3;

  a11 = 1.26;  // 1.08, 1.16
  a12 = 1.;  // 1.08, 1.16
  a13 = 1.12;  // 1.08, 1.16
  a2 = 1.1;
  a3 = 1.;

  // Special rms radii for light nucleii

  if (A < 50.)
    {
      if     (std::abs(A-1.) < 0.5)                         return 0.89*CLHEP::fermi; // p
      else if(std::abs(A-2.) < 0.5)                         return 2.13*CLHEP::fermi; // d
      else if(  // std::abs(Z-1.) < 0.5 && 
	      std::abs(A-3.) < 0.5) return 1.80*CLHEP::fermi; // t

      // else if(std::abs(Z-2.) < 0.5 && std::abs(A-3.) < 0.5) return 1.96CLHEP::fermi; // He3
      else if( // std::abs(Z-2.) < 0.5 && 
	      std::abs(A-4.) < 0.5) return 1.68*CLHEP::fermi; // He4

      else if(  // std::abs(Z-3.) < 0.5
	      std::abs(A-7.) < 0.5   )                         return 2.40*CLHEP::fermi; // Li7
      else if(  // std::abs(Z-4.) < 0.5 
	      std::abs(A-9.) < 0.5)                         return 2.51*CLHEP::fermi; // Be9

      else if( 10.  < A && A <= 16. ) r0  = a11*( 1 - (1.0/G4Pow::GetInstance()->A23(A)) )*CLHEP::fermi;   // 1.08CLHEP::fermi;
      else if( 15.  < A && A <= 20. ) r0  = a12*( 1 - (1.0/G4Pow::GetInstance()->A23(A)) )*CLHEP::fermi;
      else if( 20.  < A && A <= 30. ) r0  = a13*( 1 - (1.0/G4Pow::GetInstance()->A23(A)) )*CLHEP::fermi;
      else                            r0  = a2*CLHEP::fermi;

      R = r0*G4Pow::GetInstance()->A13(A);
    }
  else
    {
      r0 = a3*CLHEP::fermi;

      R  = r0*G4Pow::GetInstance()->powA(A, 0.27);
    }
  fNuclearRadius = R;

  return R;
}


#endif
