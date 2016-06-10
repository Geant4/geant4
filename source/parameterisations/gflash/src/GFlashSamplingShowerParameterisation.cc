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
// $Id: GFlashSamplingShowerParameterisation.cc 93079 2015-10-02 14:43:41Z gcosmo $
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ------- GFlashSamplingShowerParameterisation -------
//
// Authors: E.Barberio & Joanna Weng - 11.2005
// ------------------------------------------------------------

#include <cmath>

#include "GFlashSamplingShowerParameterisation.hh"
#include "GVFlashShowerParameterisation.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

GFlashSamplingShowerParameterisation::
GFlashSamplingShowerParameterisation(G4Material* aMat1, G4Material* aMat2,
                                     G4double dd1, G4double dd2,
                                     GFlashSamplingShowerTuning* aPar)
  : GVFlashShowerParameterisation(),
    ParAveT2(0.), ParSigLogT1(0.), ParSigLogT2(0.),
    ParSigLogA1(0.), ParSigLogA2(0.), ParRho1(0.), ParRho2(0.), ParsAveA2(0.),
    AveLogAlphah(0.), AveLogTmaxh(0.), SigmaLogAlphah(0.), SigmaLogTmaxh(0.),
    Rhoh(0.), Alphah(0.), Tmaxh(0.), Betah(0.), AveLogAlpha(0.), AveLogTmax(0.),
    SigmaLogAlpha(0.), SigmaLogTmax(0.), Rho(0.), Alpha(0.), Tmax(0.), Beta(0.)
{  
  if(!aPar) { thePar = new GFlashSamplingShowerTuning; owning = true; }
  else      { thePar = aPar; owning = false; }

  SetMaterial(aMat1,aMat2 );
  d1=dd1;
  d2=dd2;

  // Longitudinal Coefficients for a homogenious calo

  // shower max
  ParAveT1    = thePar->ParAveT1();   // ln (ln y -0.812)  
  ParAveA1    = thePar->ParAveA1();   // ln a (0.81 + (0.458 + 2.26/Z)ln y)
  ParAveA2    = thePar->ParAveA2();
  ParAveA3    = thePar->ParAveA3();
  // Sampling
  ParsAveT1   = thePar->ParsAveT1();  // T_sam = log(exp( log T_hom) + t1*Fs-1 + t2*(1-ehat));
  ParsAveT2   = thePar->ParsAveT2();
  ParsAveA1   = thePar->ParsAveA1();  
  // Variance of shower max sampling 
  ParsSigLogT1 = thePar->ParSigLogT1();    // Sigma T1 (-2.5 + 1.25 ln y)**-1 
  ParsSigLogT2 = thePar->ParSigLogT2();
  // variance of 'alpha'
  ParsSigLogA1 = thePar->ParSigLogA1();    // Sigma a (-0.82 + 0.79 ln y)**-1 
  ParsSigLogA2 = thePar->ParSigLogA2();  
  // correlation alpha%T
  ParsRho1     = thePar->ParRho1();   // Rho = 0.784 -0.023 ln y 
  ParsRho2     = thePar->ParRho2();

  // Radial Coefficients
  // r_C (tau)= z_1 +z_2 tau
  // r_t (tau)= k1 (std::exp (k3(tau -k2 ))+std::exp (k_4 (tau- k_2))))
  ParRC1 =   thePar->ParRC1();  // z_1 = 0.0251 + 0.00319 ln E
  ParRC2 =   thePar->ParRC2();
  ParRC3 =   thePar->ParRC3();  // z_2 = 0.1162 + - 0.000381 Z
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

  //additional sampling parameter
  ParsRC1= thePar->ParsRC1();
  ParsRC2= thePar->ParsRC2();
  ParsWC1= thePar->ParsWC1();
  ParsWC2=  thePar->ParsWC2(); 
  ParsRT1= thePar->ParsRT1();
  ParsRT2= thePar->ParsRT2();

  // Coeff for fluctuedted radial  profiles for a sampling media
  ParsSpotT1   = thePar->ParSpotT1();     // T_spot = T_hom =(0.698 + 0.00212)
  ParsSpotT2   = thePar->ParSpotT2();
  ParsSpotA1   = thePar->ParSpotA1();     // a_spot= a_hom (0.639 + 0.00334)
  ParsSpotA2   = thePar->ParSpotA2();    
  ParsSpotN1   = thePar->ParSpotN1();     // N_Spot 93 * ln(Z) E ** 0.876   
  ParsSpotN2   = thePar->ParSpotN2(); 
  SamplingResolution  = thePar->SamplingResolution();
  ConstantResolution  = thePar->ConstantResolution(); 
  NoiseResolution     = thePar->NoiseResolution();

  // Inits
  NSpot         = 0.00;
  AlphaNSpot    = 0.00;
  TNSpot        = 0.00;
  BetaNSpot     = 0.00;
  RadiusCore    = 0.00;
  WeightCore    = 0.00;
  RadiusTail    = 0.00;   
  ComputeZAX0EFFetc();

  G4cout << "/********************************************/ " << G4endl;
  G4cout << "  - GFlashSamplingShowerParameterisation::Constructor -  " << G4endl;
  G4cout << "/********************************************/ " << G4endl;  
}

// ------------------------------------------------------------

GFlashSamplingShowerParameterisation::~GFlashSamplingShowerParameterisation()
{
  if(owning) { delete thePar; }
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::
SetMaterial(G4Material *mat1, G4Material *mat2)
{
  G4double Es = 21*MeV;
  material1= mat1;
  Z1 = GetEffZ(material1);
  A1 = GetEffA(material1);
  density1 = material1->GetDensity();
  X01  = material1->GetRadlen();   
  Ec1      = 2.66 * std::pow((X01 * Z1 / A1),1.1); 
  Rm1 = X01*Es/Ec1;

  material2= mat2;
  Z2 = GetEffZ(material2);
  A2 = GetEffA(material2);
  density2 = material2->GetDensity();
  X02  = material2->GetRadlen();   
  Ec2      = 2.66 * std::pow((X02 * Z2 / A2),1.1); 
  Rm2 = X02*Es/Ec2; 
  // PrintMaterial(); 
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::ComputeZAX0EFFetc()
{
  G4cout << "/************ ComputeZAX0EFFetc ************/" << G4endl;
  G4cout << "  - GFlashSamplingShowerParameterisation::Material -  " << G4endl;

  G4double Es = 21*MeV; //constant

  // material and geometry parameters for a sampling calorimeter
  G4double denominator = (d1*density1 + d2*density2);
  G4double W1  = (d1*density1) / denominator;
  G4double W2  = (d2*density2) / denominator;
  Zeff   = ( W1*Z1 ) + ( W2*Z2 );    //X0*Es/Ec;
  Aeff   = ( W1*A1 ) + ( W2*A2 );
  X0eff  = ( 1./ ( ( W1 / X01) +( W2 / X02) ) ); 
  Rhoeff = ( (d1 *density1 ) + (d2 * density2 ))/G4double (d2  + d1  );
  Rmeff =  1/  ((((W1*Ec1)/ X01)   +   ((W2* Ec2)/  X02) ) / Es ) ;
  Eceff =  X0eff *((W1*Ec1)/ X01 + (W2* Ec2)/  X02 );      
  Fs =  X0eff/G4double ((d1/mm )+(d2/mm) );
  ehat = (1. / (1+ 0.007*(Z1- Z2)));

  G4cout << "W1= "  << W1 << G4endl;
  G4cout << "W2= " << W2 << G4endl;
  G4cout << "effective quantities Zeff = "<<Zeff<< G4endl;
  G4cout << "effective quantities Aeff = "<<Aeff<< G4endl;
  G4cout << "effective quantities Rhoeff = "<<Rhoeff/g *cm3<<" g/cm3"  << G4endl;
  G4cout << "effective quantities X0eff = "<<X0eff/cm <<" cm" << G4endl;

  X0eff = X0eff * Rhoeff;

  G4cout << "effective quantities X0eff = "<<X0eff/g*cm2 <<" g/cm2" << G4endl;  
  X0eff = X0eff /Rhoeff;
  G4cout << "effective quantities RMeff = "<<Rmeff/cm<<"  cm" << G4endl; 
  Rmeff = Rmeff* Rhoeff;
  G4cout << "effective quantities RMeff = "<<Rmeff/g *cm2<<" g/cm2" << G4endl;  
  Rmeff = Rmeff/ Rhoeff;
  G4cout << "effective quantities Eceff = "<<Eceff/MeV<< " MeV"<< G4endl;  
  G4cout << "effective quantities Fs = "<<Fs<<G4endl;
  G4cout << "effective quantities ehat = "<<ehat<<G4endl;
  G4cout << "/********************************************/ " <<G4endl; 
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::
GenerateLongitudinalProfile(G4double Energy)
{
  if ((material1==0) || (material2 ==0))
  {
    G4Exception("GFlashSamplingShowerParameterisation::GenerateLongitudinalProfile()",
                "InvalidSetup", FatalException, "No material initialized!");
  }  
  G4double y = Energy/Eceff;
  ComputeLongitudinalParameters(y);  
  GenerateEnergyProfile(y);
  GenerateNSpotProfile(y);
}

// ------------------------------------------------------------

void
GFlashSamplingShowerParameterisation::ComputeLongitudinalParameters(G4double y)
{
  AveLogTmaxh  = std::log(std::max(ParAveT1 +std::log(y),0.1));  //ok 
  AveLogAlphah = std::log(std::max(ParAveA1 + (ParAveA2+ParAveA3/Zeff)*std::log(y),.1)); //ok
  //hom  
  SigmaLogTmaxh  = std::min(0.5,1.00/( ParSigLogT1 + ParSigLogT2*std::log(y)) );  //ok
  SigmaLogAlphah = std::min(0.5,1.00/( ParSigLogA1 + ParSigLogA2*std::log(y)));  //ok
  Rhoh           = ParRho1+ParRho2*std::log(y);//ok
  // if sampling 
  AveLogTmax  = std::max(0.1,std::log(std::exp(AveLogTmaxh)
              + ParsAveT1/Fs + ParsAveT2*(1-ehat)));  //ok
  AveLogAlpha = std::max(0.1,std::log(std::exp(AveLogAlphah)
              + (ParsAveA1/Fs)));  //ok
  //
  SigmaLogTmax  = std::min(0.5,1.00/( ParsSigLogT1
                + ParsSigLogT2*std::log(y)) );  //ok
  SigmaLogAlpha = std::min(0.5,1.00/( ParsSigLogA1
                + ParsSigLogA2*std::log(y))); //ok
  Rho           = ParsRho1+ParsRho2*std::log(y); //ok
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::GenerateEnergyProfile(G4double /* y */)
{ 
  G4double Correlation1 = std::sqrt((1+Rho)/2);
  G4double Correlation2 = std::sqrt((1-Rho)/2);
  G4double Correlation1h = std::sqrt((1+Rhoh)/2);
  G4double Correlation2h = std::sqrt((1-Rhoh)/2);
  G4double Random1 = G4RandGauss::shoot();
  G4double Random2 = G4RandGauss::shoot();

  Tmax  = std::max(1.,std::exp( AveLogTmax  + SigmaLogTmax  *
  (Correlation1*Random1 + Correlation2*Random2) ));
  Alpha = std::max(1.1,std::exp( AveLogAlpha + SigmaLogAlpha *
  (Correlation1*Random1 - Correlation2*Random2) ));
  Beta  = (Alpha-1.00)/Tmax;
  //Parameters for Enenrgy Profile including correaltion and sigmas  
  Tmaxh  = std::exp( AveLogTmaxh  + SigmaLogTmaxh  *
  (Correlation1h*Random1 + Correlation2h*Random2) );  
  Alphah = std::exp( AveLogAlphah + SigmaLogAlphah *
  (Correlation1h*Random1 - Correlation2h*Random2) );
  Betah  = (Alphah-1.00)/Tmaxh;
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::GenerateNSpotProfile(const G4double y)
{
  TNSpot     = Tmaxh *  (ParsSpotT1+ParsSpotT2*Zeff); //ok.
  TNSpot     = std::max(0.5,Tmaxh * (ParsSpotT1+ParsSpotT2*Zeff));
  AlphaNSpot = Alphah * (ParsSpotA1+ParsSpotA2*Zeff);     
  BetaNSpot  = (AlphaNSpot-1.00)/TNSpot;           // ok 
  NSpot      = ParsSpotN1 /SamplingResolution * std::pow(y*Eceff/GeV,ParsSpotN2 );
}

// ------------------------------------------------------------

G4double
GFlashSamplingShowerParameterisation::
ApplySampling(const G4double DEne, const G4double )
{
  G4double DEneFluctuated = DEne;
  G4double Resolution     = std::pow(SamplingResolution,2);

  //       +pow(NoiseResolution,2)/  //@@@@@@@@ FIXME 
  //                         Energy*(1.*MeV)+
  //                         pow(ConstantResolution,2)*
  //                          Energy/(1.*MeV);

  if(Resolution >0.0 && DEne > 0.00)
  {
    G4float x1=DEne/Resolution;
    G4float x2 = G4RandGamma::shoot(x1, 1.0)*Resolution;     
    DEneFluctuated=x2;
  }
  return DEneFluctuated;
}

// ------------------------------------------------------------

G4double GFlashSamplingShowerParameterisation::
IntegrateEneLongitudinal(G4double LongitudinalStep)
{
  G4double LongitudinalStepInX0 = LongitudinalStep / X0eff;
  G4float x1= Betah*LongitudinalStepInX0;
  G4float x2= Alphah;
  float x3 =  gam(x1,x2);
  G4double DEne=x3;
  return DEne;
}

// ------------------------------------------------------------

G4double GFlashSamplingShowerParameterisation::
IntegrateNspLongitudinal(G4double LongitudinalStep)
{
  G4double LongitudinalStepInX0 = LongitudinalStep / X0eff; 
  G4float x1 = BetaNSpot*LongitudinalStepInX0;
  G4float x2 = AlphaNSpot;
  G4float x3 =  gam(x1,x2);
  G4double DNsp = x3;
  return DNsp;
}

// ------------------------------------------------------------

G4double GFlashSamplingShowerParameterisation::
GenerateRadius(G4int ispot, G4double Energy, G4double LongitudinalPosition)
{
  if(ispot < 1) 
  {
    // Determine lateral parameters in the middle of the step.
    // They depend on energy & position along step
    //
    G4double Tau = ComputeTau(LongitudinalPosition);
    ComputeRadialParameters(Energy,Tau);  
  }
  
  G4double Radius;
  G4double Random1 = G4UniformRand();
  G4double Random2 = G4UniformRand(); 
  if(Random1  <WeightCore) //WeightCore = p < w_i  
  {
    Radius = Rmeff * RadiusCore * std::sqrt( Random2/(1. - Random2) );
  }
  else
  {
    Radius = Rmeff * RadiusTail * std::sqrt( Random2/(1. - Random2) );
  }   
  Radius =  std::min(Radius,DBL_MAX);
  return Radius;
}

// ------------------------------------------------------------

G4double
GFlashSamplingShowerParameterisation::
ComputeTau(G4double LongitudinalPosition)
{
  G4double tau = LongitudinalPosition / Tmax/ X0eff     //<t> = T* a /(a - 1) 
                 * (Alpha-1.00) /Alpha
                 * std::exp(AveLogAlpha)/(std::exp(AveLogAlpha)-1.);  //ok 
  return tau;
}

// ------------------------------------------------------------

void GFlashSamplingShowerParameterisation::
ComputeRadialParameters(G4double Energy, G4double Tau)
{
  G4double z1 = ParRC1 + ParRC2* std::log(Energy/GeV);         //ok
  G4double z2 = ParRC3+ParRC4*Zeff;                            //ok
  RadiusCore  =  z1 + z2 * Tau;                                //ok 
  G4double p1 = ParWC1+ParWC2*Zeff;                            //ok
  G4double p2 = ParWC3+ParWC4*Zeff;                            //ok
  G4double p3 = ParWC5+ParWC6*std::log(Energy/GeV);            //ok
  WeightCore   =  p1 * std::exp( (p2-Tau)/p3-  std::exp( (p2-Tau) /p3) ); //ok
  
  G4double k1 = ParRT1+ParRT2*Zeff;                // ok
  G4double k2 = ParRT3;                            // ok
  G4double k3 = ParRT4;                            // ok
  G4double k4 = ParRT5+ParRT6* std::log(Energy/GeV);    // ok
  
  RadiusTail   = k1*(std::exp(k3*(Tau-k2))
               + std::exp(k4*(Tau-k2)) );            //ok

  // sampling calorimeter  

  RadiusCore   = RadiusCore + ParsRC1*(1-ehat) + ParsRC2/Fs*std::exp(-Tau); //ok
  WeightCore   = WeightCore + (1-ehat)
                            * (ParsWC1+ParsWC2/Fs * std::exp(-std::pow((Tau-1.),2))); //ok
  RadiusTail   = RadiusTail + (1-ehat)* ParsRT1+ ParsRT2/Fs *std::exp(-Tau);     //ok  
}

// ------------------------------------------------------------

G4double GFlashSamplingShowerParameterisation::
GenerateExponential(const G4double /* Energy */ )
{
  G4double ParExp1 =  9./7.*X0eff;
  G4double random  = -ParExp1*G4RandExponential::shoot() ;
  return random;
}
