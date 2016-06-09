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
// $Id: G4DiffuseElastic.hh,v 1.7 2007/06/12 14:46:26 grichine Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// G4 Model: optical elastic scattering with 4-momentum balance 
//
// Class Description
// Final state production model for hadron nuclear elastic scattering; 
// Class Description - End
//
//
// 24.05.07 V. Grichine first implementation for proton elastic scattering


#ifndef G4DiffuseElastic_h
#define G4DiffuseElastic_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

using namespace std;

class G4ParticleDefinition;

class G4DiffuseElastic : public G4HadronicInteraction
{
public:

  G4DiffuseElastic();

  virtual ~G4DiffuseElastic();
 
  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
				  G4Nucleus & targetNucleus);


  void SetPlabLowLimit(G4double value);

  void SetHEModelLowLimit(G4double value);

  void SetQModelLowLimit(G4double value);

  void SetLowestEnergyLimit(G4double value);

  void SetRecoilKinEnergyLimit(G4double value);

  G4double SampleT(const G4ParticleDefinition* aParticle, 
                         G4double p, G4double A);

  G4double SampleThetaCMS(const G4ParticleDefinition* aParticle, G4double p, G4double A);

  G4double SampleThetaLab(const G4HadProjectile* aParticle, 
                                G4double tmass, G4double A);

  G4double GetDiffuseElasticXsc( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A         );
  G4double IntegralElasticProb( const G4ParticleDefinition* particle, 
                                 G4double theta, 
			         G4double momentum, 
				 G4double A            );
  
  G4double BesselJzero(G4double z);
  G4double BesselJone(G4double z);
  G4double DampFactor(G4double z);
  G4double BesselOneByArg(G4double z);

  G4double GetDiffElasticProb(G4double theta);
  G4double GetIntegrandFunction(G4double theta);


  G4double GetNuclearRadius(){return fNuclearRadius;};

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

  const G4ParticleDefinition* fParticle;
  G4double fWaveVector;
  G4double fAtomicWeight;
  G4double fNuclearRadius;

};

inline void G4DiffuseElastic::SetRecoilKinEnergyLimit(G4double value)
{
  lowEnergyRecoilLimit = value;
}

inline void G4DiffuseElastic::SetPlabLowLimit(G4double value)
{
  plabLowLimit = value;
}

inline void G4DiffuseElastic::SetHEModelLowLimit(G4double value)
{
  lowEnergyLimitHE = value;
}

inline void G4DiffuseElastic::SetQModelLowLimit(G4double value)
{
  lowEnergyLimitQ = value;
}

inline void G4DiffuseElastic::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}


/////////////////////////////////////////////////////////////
//
// Bessel J0 function based on rational approximation from 
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141 

inline G4double G4DiffuseElastic::BesselJzero(G4double value)
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

inline G4double G4DiffuseElastic::BesselJone(G4double value)
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
// damp factor in diffraction x*pi/sh(x*pi)

inline G4double G4DiffuseElastic::DampFactor(G4double x)
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

inline G4double G4DiffuseElastic::BesselOneByArg(G4double x)
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

#endif
