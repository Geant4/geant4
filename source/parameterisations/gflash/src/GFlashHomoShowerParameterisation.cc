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
// $Id: GFlashHomoShowerParameterisation.cc 69579 2013-05-08 13:53:57Z gcosmo $
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ------- GFlashHomoShowerParameterisation -------
//
// Authors: E.Barberio & Joanna Weng - 9.11.2004
// ------------------------------------------------------------

#include <cmath>

#include "GFlashHomoShowerParameterisation.hh"
#include "GVFlashShowerParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

GFlashHomoShowerParameterisation::
GFlashHomoShowerParameterisation(G4Material * aMat,
                                 GVFlashHomoShowerTuning * aPar)
  : GVFlashShowerParameterisation(),
    ConstantResolution(0.), NoiseResolution(0.), SamplingResolution(0.),
    AveLogAlphah(0.), AveLogTmaxh(0.), SigmaLogAlphah(0.), SigmaLogTmaxh(0.),
    Rhoh(0.), Alphah(0.), Tmaxh(0.), Betah(0.)

{  
  if(!aPar) { thePar = new GVFlashHomoShowerTuning; owning = true; }
  else      { thePar = aPar; owning = false; }

  SetMaterial(aMat);
  PrintMaterial(aMat);

  /********************************************/
  /* Homo Calorimeter                         */
  /********************************************/ 
  // Longitudinal Coefficients for a homogenious calo
  // shower max
  //
  ParAveT1    = thePar->ParAveT1();   // ln (ln y -0.812)  
  ParAveA1    = thePar->ParAveA1();   // ln a (0.81 + (0.458 + 2.26/Z)ln y)
  ParAveA2    = thePar->ParAveA2();
  ParAveA3    = thePar->ParAveA3();

  // Variance of shower max
  ParSigLogT1 = thePar->ParSigLogT1();     // Sigma T1 (-1.4 + 1.26 ln y)**-1 
  ParSigLogT2 = thePar->ParSigLogT2();

  // variance of 'alpha'
  //
  ParSigLogA1 = thePar->ParSigLogA1();    // Sigma a (-0.58 + 0.86 ln y)**-1 
  ParSigLogA2 = thePar->ParSigLogA2();

  // correlation alpha%T
  //
  ParRho1     = thePar->ParRho1();   // Rho = 0.705 -0.023 ln y 
  ParRho2     = thePar->ParRho2();

  // Radial Coefficients
  // r_C (tau)= z_1 +z_2 tau
  // r_t (tau)= k1 (std::exp (k3(tau -k2 ))+std::exp (k_4 (tau- k_2))))
  //
  ParRC1 =   thePar->ParRC1();   // z_1 = 0.0251 + 0.00319 ln E
  ParRC2 =   thePar->ParRC2();

  ParRC3 =   thePar->ParRC3();   // z_2 = 0.1162 + - 0.000381 Z
  ParRC4 =   thePar->ParRC4();

  ParWC1 = thePar->ParWC1();
  ParWC2 = thePar->ParWC2();
  ParWC3 = thePar->ParWC3();
  ParWC4 = thePar->ParWC4();
  ParWC5 = thePar->ParWC5(); 
  ParWC6 = thePar->ParWC6();

  ParRT1 = thePar->ParRT1();
  ParRT2 = thePar->ParRT2();
  ParRT3 = thePar->ParRT3();
  ParRT4 = thePar->ParRT4(); 
  ParRT5 = thePar->ParRT5();
  ParRT6 = thePar->ParRT6();

  // Coeff for fluctueted radial  profiles for a uniform media
  //
  ParSpotT1   = thePar->ParSpotT1();     // T_spot = T_hom =(0.698 + 0.00212)
  ParSpotT2   = thePar->ParSpotT2();

  ParSpotA1   = thePar->ParSpotA1();     // a_spot= a_hom (0.639 + 0.00334)
  ParSpotA2   = thePar->ParSpotA2();  

  ParSpotN1   = thePar->ParSpotN1();     // N_Spot 93 * ln(Z) E ** 0.876   
  ParSpotN2   = thePar->ParSpotN2(); 

  // Inits

  NSpot         = 0.00;
  AlphaNSpot    = 0.00;
  TNSpot        = 0.00;
  BetaNSpot     = 0.00;

  RadiusCore    = 0.00;
  WeightCore    = 0.00;
  RadiusTail    = 0.00; 

  G4cout << "/********************************************/ " << G4endl;
  G4cout << "  - GFlashHomoShowerParameterisation::Constructor -  " << G4endl;  
  G4cout << "/********************************************/ " << G4endl;
}

void GFlashHomoShowerParameterisation::SetMaterial(G4Material *mat)
{
  material= mat;
  Z = GetEffZ(material);
  A = GetEffA(material);
  density = material->GetDensity()/(g/cm3);
  X0  = material->GetRadlen(); 
  Ec      = 2.66 * std::pow((X0 * Z / A),1.1); 
  G4double Es = 21*MeV;
  Rm = X0*Es/Ec;
  // PrintMaterial(); 
}

GFlashHomoShowerParameterisation::~GFlashHomoShowerParameterisation()
{
  if(owning) { delete thePar; }
}

void GFlashHomoShowerParameterisation::
GenerateLongitudinalProfile(G4double Energy)
{
  if (material==0) 
  {
    G4Exception("GFlashHomoShowerParameterisation::GenerateLongitudinalProfile()",
                "InvalidSetup", FatalException, "No material initialized!");
  }
  
  G4double y = Energy/Ec;
  ComputeLongitudinalParameters(y);  
  GenerateEnergyProfile(y);
  GenerateNSpotProfile(y);
}

void
GFlashHomoShowerParameterisation::ComputeLongitudinalParameters(G4double y)
{
  AveLogTmaxh  = std::log(ParAveT1 + std::log(y));
    //ok  <ln T hom>
  AveLogAlphah = std::log(ParAveA1 + (ParAveA2+ParAveA3/Z)*std::log(y));
    //ok  <ln alpha hom> 

  SigmaLogTmaxh  = 1.00/( ParSigLogT1 + ParSigLogT2*std::log(y)) ;
    //ok sigma (ln T hom)
  SigmaLogAlphah = 1.00/( ParSigLogA1 + ParSigLogA2*std::log(y));
    //ok sigma (ln alpha hom)
  Rhoh           = ParRho1+ParRho2*std::log(y);        //ok             
}

void GFlashHomoShowerParameterisation::GenerateEnergyProfile(G4double /* y */)
{ 
  G4double Correlation1h = std::sqrt((1+Rhoh)/2);
  G4double Correlation2h = std::sqrt((1-Rhoh)/2);

  G4double Random1 = G4RandGauss::shoot();
  G4double Random2 = G4RandGauss::shoot();

  // Parameters for Enenrgy Profile including correaltion and sigmas 
 
  Tmaxh  = std::exp( AveLogTmaxh  + SigmaLogTmaxh  *
           (Correlation1h*Random1 + Correlation2h*Random2) );
  Alphah = std::exp( AveLogAlphah + SigmaLogAlphah *
           (Correlation1h*Random1 - Correlation2h*Random2) );
  Betah  = (Alphah-1.00)/Tmaxh;
}

void GFlashHomoShowerParameterisation::GenerateNSpotProfile(const G4double y)
{
  TNSpot     = Tmaxh *  (ParSpotT1+ParSpotT2*Z);   // ok
  AlphaNSpot = Alphah * (ParSpotA1+ParSpotA2*Z);   
  BetaNSpot  = (AlphaNSpot-1.00)/TNSpot;           // ok 
  NSpot      = ParSpotN1 * std::log(Z)*std::pow((y*Ec)/GeV,ParSpotN2 ); // ok
}

G4double GFlashHomoShowerParameterisation::
IntegrateEneLongitudinal(G4double LongitudinalStep)
{
  G4double LongitudinalStepInX0 = LongitudinalStep / X0;
  G4float x1= Betah*LongitudinalStepInX0;
  G4float x2= Alphah;
  float x3 =  gam(x1,x2);
  G4double DEne=x3;
  return DEne;
}

G4double GFlashHomoShowerParameterisation::
IntegrateNspLongitudinal(G4double LongitudinalStep)
{
  G4double LongitudinalStepInX0 = LongitudinalStep / X0; 
  G4float x1 = BetaNSpot*LongitudinalStepInX0;
  G4float x2 = AlphaNSpot;
  G4float x3 =  gam(x1,x2);
  G4double DNsp = x3;
  return DNsp;
}


G4double GFlashHomoShowerParameterisation::
GenerateRadius(G4int ispot, G4double Energy, G4double LongitudinalPosition)
{
  if(ispot < 1) 
  {
    // Determine lateral parameters in the middle of the step.
    // They depend on energy & position along step.
    //
    G4double Tau = ComputeTau(LongitudinalPosition);
    ComputeRadialParameters(Energy,Tau);  
  }

  G4double Radius;
  G4double Random1 = G4UniformRand();
  G4double Random2 = G4UniformRand(); 

  if(Random1  <WeightCore) //WeightCore = p < w_i  
  {
    Radius = Rm * RadiusCore * std::sqrt( Random2/(1. - Random2) );
  }
  else
  {
    Radius = Rm * RadiusTail * std::sqrt( Random2/(1. - Random2) );
  }   
  Radius =  std::min(Radius,DBL_MAX);
  return Radius;
}

G4double GFlashHomoShowerParameterisation::
ComputeTau(G4double LongitudinalPosition)
{
  G4double tau = LongitudinalPosition / Tmaxh / X0     //<t> = T* a /(a - 1) 
  * (Alphah-1.00) /Alphah * 
  std::exp(AveLogAlphah)/(std::exp(AveLogAlphah)-1.);  //ok 
  return tau;
}

void GFlashHomoShowerParameterisation::
ComputeRadialParameters(G4double Energy, G4double Tau)
{
  G4double z1 = ParRC1 + ParRC2* std::log(Energy/GeV)  ;       //ok
  G4double z2 = ParRC3+ParRC4*Z ;                              //ok
  RadiusCore   =  z1 + z2 * Tau  ;                             //ok 

  G4double p1 = ParWC1+ParWC2*Z;                               //ok
  G4double p2 = ParWC3+ParWC4*Z;                               //ok
  G4double p3 = ParWC5+ParWC6*std::log(Energy/GeV);            //ok

  WeightCore   =  p1 * std::exp( (p2-Tau)/p3 - std::exp( (p2-Tau) /p3) ); //ok

  G4double k1 = ParRT1+ParRT2*Z;                   // ok
  G4double k2 = ParRT3;                            // ok
  G4double k3 = ParRT4;                            // ok
  G4double k4 = ParRT5+ParRT6* std::log(Energy/GeV);    // ok

  RadiusTail   = k1*(std::exp(k3*(Tau-k2))  +
  std::exp(k4*(Tau-k2)) );            //ok
}

G4double GFlashHomoShowerParameterisation::
GenerateExponential(const G4double /* Energy */ )
{
  G4double ParExp1 =  9./7.*X0;
  G4double random  = -ParExp1*G4RandExponential::shoot() ;
  return random;
}
