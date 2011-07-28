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
// $Id: G4AntiNuclElastic.cc  - A.Galoyan 02.05.2011
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Geant4 Header : G4AntiNuclElastic
//
// 

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiAlpha.hh"
#include "G4AntiTriton.hh"
#include "G4AntiHe3.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4Pow.hh"
#include "G4AntiNuclElastic.hh"

#include "G4NucleiProperties.hh"       // 8 May 2011

G4AntiNuclElastic::G4AntiNuclElastic() 
  : G4HadronElastic("AntiAElastic")
{

  SetMinEnergy( 0.1*GeV );
  SetMaxEnergy( 1.*TeV );


  theAProton       = G4AntiProton::AntiProton();
  theANeutron      = G4AntiNeutron::AntiNeutron();
  theADeuteron     = G4AntiDeuteron::AntiDeuteron();
  theATriton       = G4AntiTriton::AntiTriton();
  theAAlpha        = G4AntiAlpha::AntiAlpha();
  theAHe3          = G4AntiHe3::AntiHe3();

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();


  cs = new G4ComponentAntiNuclNuclearXS();
  fParticle = 0;
  fWaveVector = 0.;
  fBeta = 0.;
  fZommerfeld = 0.;
  fAm = 0.;
  fTetaCMS = 0.;    
  fRa = 0.;
  fRef = 0.;
  fceff = 0.;
  fptot = 0.;
  fTmax = 0.;
  fThetaLab = 0.;
}

/////////////////////////////////////////////////////////////////////////
G4AntiNuclElastic::~G4AntiNuclElastic()
{
  delete cs;  // 8 May 2011
}

////////////////////////////////////////////////////////////////////////
// sample momentum transfer in the CMS system 
G4double G4AntiNuclElastic::SampleInvariantT(const G4ParticleDefinition* particle, 
				   G4double Plab,  G4int Z, G4int A)
{
   G4double T;
   G4double Mproj = particle->GetPDGMass();
   G4LorentzVector Pproj(0.,0.,Plab,std::sqrt(Plab*Plab+Mproj*Mproj));
   G4double ctet1 = GetcosTeta1(Plab, A); 

   G4double energy=Pproj.e()-Mproj;   

  const G4ParticleDefinition* theParticle = particle;

  G4ParticleDefinition * theDef = 0;

  if(Z == 1 && A == 1)       theDef = theProton;
  else if (Z == 1 && A == 2) theDef = theDeuteron;
  else if (Z == 1 && A == 3) theDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theDef = G4He3::He3();
  else if (Z == 2 && A == 4) theDef = theAlpha;
/*                                                            // 8 May 2011
  else theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);

  G4double TargMass = theDef->GetPDGMass(); 
*/                                                            // 8 May 2011

  G4double TargMass =G4NucleiProperties::GetNuclearMass(A,Z); // 8 May 2011

  //transform to CMS

   G4LorentzVector lv(0.0,0.0,0.0,TargMass);   
   lv += Pproj;
  G4double S = lv.mag2()/GeV/GeV;
  G4ThreeVector bst = lv.boostVector();
  Pproj.boost(-bst);

  G4ThreeVector p1 = Pproj.vect();
  G4double ptot    = p1.mag();

  fbst = bst;
  fptot= ptot;
  fTmax = 4.0*ptot*ptot;  
  
  G4double  Z1 = particle->GetPDGCharge();
  G4double  Z2 = Z;  
  
   G4double beta = CalculateParticleBeta(particle, ptot); 
   G4double n  = CalculateZommerfeld( beta,  Z1,  Z2 );
   G4double Am = CalculateAm( ptot,  n,  Z2 );
   fWaveVector = ptot;     //   /hbarc; 
    
   G4LorentzVector Fproj(0.,0.,0.,0.);  
 
 
G4double  XsCoulomb = sqr(n/fWaveVector)*pi*(1+ctet1)/(1.+Am)/(1.+2.*Am-ctet1); 
XsCoulomb=XsCoulomb*0.38938e+6;
//XsCoulomb=0;

// 8 May 2011
//G4double XsElastHad =cs->ComputeElasticCrossSection(particle, energy, Z, (G4double)A);
//G4double XstotalHad =cs->ComputeTotalCrossSection(particle, energy, Z, (G4double)A);

G4double XsElastHad =cs->GetElasticElementCrossSection(particle, energy, Z, (G4double)A);
G4double XstotalHad =cs->GetTotalElementCrossSection(particle, energy, Z, (G4double)A);
// 8 May 2011

//XsElastHad=0.;

 XsElastHad/=millibarn; XstotalHad/=millibarn;
   G4double CoulombProb =  XsCoulomb/(XsCoulomb+XsElastHad);

//G4cout<<" XsCol"<< XsCoulomb <<" Elast" << XsElastHad<<G4endl;

  if(G4UniformRand() < CoulombProb)
  {  // Simulation of Coulomb scattering

   G4double phi = twopi * G4UniformRand();
   G4double Ksi =  G4UniformRand();

   G4double par1 = 2.*(1.+Am)/(1.+ctet1);

// ////sample ThetaCMS in Coulomb part

   G4double cosThetaCMS = (par1*ctet1- Ksi*(1.+2.*Am))/(par1-Ksi);

   G4double PtZ=ptot*cosThetaCMS;
   Fproj.setPz(PtZ);
   G4double PtProjCMS = ptot*std::sqrt(1.0 - cosThetaCMS*cosThetaCMS);
   G4double PtX= PtProjCMS * std::cos(phi);
   G4double PtY= PtProjCMS * std::sin(phi);
   Fproj.setPx(PtX);     
   Fproj.setPy(PtY);
   Fproj.setE(std::sqrt(PtX*PtX+PtY*PtY+PtZ*PtZ+Mproj*Mproj));    
   T =  -(Pproj-Fproj).mag2();      
  } else

 {  
///////Simulation of strong interaction scattering////////////////////////////

   G4double Qmax = 2.*ptot*197.33;    // in fm^-1
   G4double Amag = 70*70;    //  A1 in Magora funct:A1*exp(-q*A2)
   G4double SlopeMag = 2.*3.0; // A2 in Magora funct:A1*exp(-q*A2)

  G4double sig_pbarp= cs->GetAntiHadronNucleonTotCrSc(particle,energy); 
   
// fRa = 1.113*std::pow(A,0.333333) - 0.227/std::pow(A,0.333333); // 8 May 2011
   fRa = 1.113*G4Pow::GetInstance()->Z13(A) -                     // 8 May 2011
         0.227/G4Pow::GetInstance()->Z13(A);                      // 8 May 2011
   if(A == 3) fRa=1.81;
   if(A == 4) fRa=1.37;
 
       if((A>=12.) && (A<27) ) fRa=fRa*0.85;
       if((A>=27.) && (A<48) ) fRa=fRa*0.90;
       if((A>=48.) && (A<65) ) fRa=fRa*0.95;

 G4double Ref2 = 0;
 G4double ceff2 =0;
 // G4double A2 = 0;

 if  ((theParticle == theAProton) || (theParticle == theANeutron))  
{
  if(theDef == theProton)
  { 
   G4double Mp2=sqr(theDef->GetPDGMass()/GeV ); 
// Ref2=0.8/pow((S-4.*Mp2),0.25) + 0.55;           // 8 May 2011
   Ref2=0.8/std::sqrt(std::sqrt(S-4.*Mp2)) + 0.55; // 8 May 2011
   Ref2=Ref2*Ref2;
   ceff2 = 0.035/(sqr(S-4.3)+0.4) + 0.085 * std::log(S) ;
   ceff2 = ceff2*ceff2;
   //A2 = (4.0/std::sqrt(S-4.*Mp2) + 0.3) / 256.8/pi;

SlopeMag = 0.5;  // Uzhi
Amag=1.;         // Uzhi
  }
 if(Z>2) 
 {  Ref2 = fRa*fRa +2.48*0.01*sig_pbarp*fRa - 2.23e-6*sig_pbarp*sig_pbarp*fRa*fRa; 
    ceff2 = 0.16+3.3e-4*sig_pbarp+0.35*exp(-0.03*sig_pbarp);
 }
 if( (Z==2)&&(A==4) )
 {  Ref2 = fRa*fRa -0.46 +0.03*sig_pbarp - 2.98e-6*sig_pbarp*sig_pbarp;
    ceff2= 0.078 + 6.657e-4*sig_pbarp + 0.3359*exp(-0.03*sig_pbarp);
 }
  if( (Z==1)&&(A==3) )
 {  Ref2 = fRa*fRa - 1.36 + 0.025 * sig_pbarp - 3.69e-7 * sig_pbarp*sig_pbarp;
    ceff2 = 0.149 + 7.091e-04*sig_pbarp + 0.3743*exp(-0.03*sig_pbarp);
 }
  if( (Z==2)&&(A==3) )
 {  Ref2 = fRa*fRa - 1.36 + 0.025 * sig_pbarp - 3.69e-7 * sig_pbarp*sig_pbarp;
    ceff2 = 0.149 + 7.091e-04*sig_pbarp + 0.3743*exp(-0.03*sig_pbarp);
 }  
  if( (Z==1)&&(A==2) )
 {
    Ref2 = fRa*fRa - 0.28 + 0.019 * sig_pbarp + 2.06e-6 * sig_pbarp*sig_pbarp;  
    ceff2 = 0.297 + 7.853e-04*sig_pbarp + 0.2899*exp(-0.03*sig_pbarp);
 }
 }

 if (theParticle == theADeuteron)
{
 sig_pbarp= cs->GetAntiHadronNucleonTotCrSc(particle,energy/2.); 
  Ref2 = XstotalHad/10./2./pi ;
  if(Z>2)
  {
   ceff2 = 0.38 + 2.0e-4 *sig_pbarp + 0.5 * exp(-0.03*sig_pbarp);
  }
  if(theDef == theProton)
  { 
    ceff2 = 0.297 + 7.853e-04*sig_pbarp + 0.2899*exp(-0.03*sig_pbarp);
  }
  if(theDef == theDeuteron)
  { 
   ceff2 = 0.65 + 3.0e-4*sig_pbarp + 0.55 * exp(-0.03*sig_pbarp);
  }
  if( (theDef == G4Triton::Triton()) || (theDef == G4He3::He3() ) )
  {
   ceff2 = 0.57 + 2.5e-4*sig_pbarp + 0.65 * exp(-0.02*sig_pbarp);
  }
  if(theDef == theAlpha)
  {
   ceff2 = 0.40 + 3.5e-4 *sig_pbarp + 0.45 * exp(-0.02*sig_pbarp);
  }
}

   if( (theParticle ==theAHe3) || (theParticle ==theATriton) )
 {
   sig_pbarp = cs->GetAntiHadronNucleonTotCrSc(particle,energy/3.);
   Ref2 = XstotalHad/10./2./pi ;
  if(Z>2)
  {
   ceff2 = 0.26 + 2.2e-4*sig_pbarp + 0.33*exp(-0.03*sig_pbarp);
  }  
    if(theDef == theProton)
  {
   ceff2 = 0.149 + 7.091e-04*sig_pbarp + 0.3743*exp(-0.03*sig_pbarp);         
  }
   if(theDef == theDeuteron)
  {
   ceff2 = 0.57 + 2.5e-4*sig_pbarp + 0.65 * exp(-0.02*sig_pbarp);
  }
 if( (theDef == G4Triton::Triton()) || (theDef == G4He3::He3() ) )
  {
    ceff2 = 0.39 + 2.7e-4*sig_pbarp + 0.7 * exp(-0.02*sig_pbarp);
  }
  if(theDef == theAlpha)
  {
   ceff2 = 0.24 + 3.5e-4*sig_pbarp + 0.75 * exp(-0.03*sig_pbarp);
  }
 }


 if (theParticle == theAAlpha)
 {
  sig_pbarp = cs->GetAntiHadronNucleonTotCrSc(particle,energy/3.);
   Ref2 = XstotalHad/10./2./pi ;
  if(Z>2)
  {
   ceff2 = 0.22 + 2.0e-4*sig_pbarp + 0.2 * exp(-0.03*sig_pbarp);
  }
    if(theDef == theProton)
  {
    ceff2= 0.078 + 6.657e-4*sig_pbarp + 0.3359*exp(-0.03*sig_pbarp);   
  }
     if(theDef == theDeuteron)  
  {
   ceff2 = 0.40 + 3.5e-4 *sig_pbarp + 0.45 * exp(-0.02*sig_pbarp);
  }   
  if( (theDef == G4Triton::Triton()) || (theDef == G4He3::He3() ) )
  {
    ceff2 = 0.24 + 3.5e-4*sig_pbarp + 0.75 * exp(-0.03*sig_pbarp);   
  }
 if(theDef == theAlpha)   
  {
   ceff2 = 0.17 + 3.5e-4*sig_pbarp + 0.45 * exp(-0.03*sig_pbarp);
  }
}

    fRef=std::sqrt(Ref2);
    fceff = std::sqrt(ceff2);     
//    G4cout<<" Ref  "<<fRef<<" c_eff "<<fceff<<G4endl; 

  
   G4double Q = 0.0 ;
  G4double BracFunct;
 do 
 {
  Q = -std::log(1.-(1.- std::exp(-SlopeMag * Qmax))* G4UniformRand() )/SlopeMag;
  G4double x = fRef * Q;

 BracFunct = sqr(BesselOneByArg(x)*DampFactor(pi*fceff*Q))/
 (Amag*std::exp(-SlopeMag*Q));

//G4cout<<"Amag Bess "<<Amag<<" "<<BesselOneByArg(0.)*DampFactor(pi*fceff*0.)<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

// BracFunct = BracFunct * Q * std::pow(fRef,4.); // 8 May 2011
   BracFunct = BracFunct * Q * sqr(sqr(fRef));   // 8 May 2011
  } 
   while (G4UniformRand()>BracFunct);
    T= sqr(Q); 
    T*=3.893913e+4;                // fm -> MeV^2
      }

 G4double cosTet=1.0-T/(2.*ptot*ptot);
 
 fTetaCMS=std::acos(cosTet);

 return T;
}

/////////////////////////////////////////////////////////////////////
//  Sample of Theta in CMS
 G4double G4AntiNuclElastic::SampleThetaCMS(const G4ParticleDefinition* p, G4double plab,
                                                                         G4int Z, G4int A)
{ 
  G4double T;
  T =  SampleInvariantT( p, plab,  Z,  A);

   // NaN finder
 if(!(T < 0.0 || T >= 0.0))
  {
    if (verboseLevel > 0)
    {
      G4cout << "G4DiffuseElastic:WARNING: A = " << A
             << " mom(GeV)= " << plab/GeV
             << " S-wave will be sampled"
             << G4endl;
    }
    T = G4UniformRand()*fTmax;
 
  }

 G4double cosTet=1.0-T/(2.*fptot*fptot);
 fTetaCMS=std::acos(cosTet); 
 return fTetaCMS;
}  


/////////////////////////////////////////////////////////////////////
//  Sample of Theta in Lab System
 G4double G4AntiNuclElastic::SampleThetaLab(const G4ParticleDefinition* p, G4double plab,
                                                                         G4int Z, G4int A)
{ 

  G4double T; 
  T = SampleInvariantT( p, plab,  Z,  A);

 // NaN finder
 if(!(T < 0.0 || T >= 0.0))
  {
    if (verboseLevel > 0)               
    {
      G4cout << "G4DiffuseElastic:WARNING: A = " << A
             << " mom(GeV)= " << plab/GeV
             << " S-wave will be sampled"
             << G4endl;
    }
    T = G4UniformRand()*fTmax;
  }

    G4double phi  = G4UniformRand()*twopi;
    G4double cost = 1. - 2.0*T/fTmax;
    G4double sint;
  if( cost >= 1.0 )
  {
    cost = 1.0;
    sint = 0.0;
  }
  else if( cost <= -1.0)
  {
    cost = -1.0;
    sint =  0.0;
  }
 else
  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }

 G4double m1 = p->GetPDGMass();
 G4ThreeVector v(sint*std::cos(phi),sint*std::sin(phi),cost);
 v *= fptot;
 G4LorentzVector nlv(v.x(),v.y(),v.z(),std::sqrt(fptot*fptot + m1*m1));
   
  nlv.boost(fbst);
   
  G4ThreeVector np = nlv.vect();
  G4double theta = np.theta();
  fThetaLab = theta; 

  return theta;
}

////////////////////////////////////////////////////////////////////
//   Calculation of Damp factor
 G4double G4AntiNuclElastic::DampFactor(G4double x)
{
   G4double df;
   G4double  f3 = 6.; // first factorials

 if( std::fabs(x) < 0.01 )
  { 
      df=1./(1.+x*x/f3); 
 }
  else
  {
    df = x/std::sinh(x); 
  }
  return df;
}


/////////////////////////////////////////////////////////////////////////////////
//  Calculation of particle velocity Beta

 G4double G4AntiNuclElastic::CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    )
{
  G4double mass = particle->GetPDGMass();
  G4double a    = momentum/mass;
  fBeta         = a/std::sqrt(1+a*a);

  return fBeta; 
}


///////////////////////////////////////////////////////////////////////////////////
//   Calculation of parameter Zommerfeld

 G4double G4AntiNuclElastic::CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 )
{
  fZommerfeld = fine_structure_const*Z1*Z2/beta;

  return fZommerfeld; 
}

////////////////////////////////////////////////////////////////////////////////////
//  
G4double G4AntiNuclElastic::CalculateAm( G4double momentum, G4double n, G4double Z)
{
  G4double k   = momentum/hbarc;
  G4double ch  = 1.13 + 3.76*n*n;
//G4double zn  = 1.77*k*std::pow(Z,-1./3.)*Bohr_radius;           // 8 May 2011
  G4double zn  = 1.77*k/G4Pow::GetInstance()->A13(Z)*Bohr_radius; // 8 May 2011
  G4double zn2 = zn*zn;
  fAm          = ch/zn2;

  return fAm;
}

//////////////////////////////////////////////////////////////////////////////
// Bessel J1 function based on rational approximation from
// J.F. Hart, Computer Approximations, New York, Willey 1968, p. 141
    
 G4double G4AntiNuclElastic::BesselJone(G4double value)
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

////////////////////////////////////////////////////////////////////////////////
// return J1(x)/x with special case for small x 
G4double G4AntiNuclElastic::BesselOneByArg(G4double x)
{
  G4double x2, result;
 
  if( std::fabs(x) < 0.01 )
  {
   x  *= 0.5;
   x2 = x*x;
   result = (2.- x2 + x2*x2/6.)/4.;
  }
  else
  {
    result = BesselJone(x)/x;
  }
  return result;
} 

/////////////////////////////////////////////////////////////////////////////////
// return  angle from which Coulomb scattering is calculated 
G4double G4AntiNuclElastic::GetcosTeta1(G4double plab, G4int A)
{

// G4double p0 =G4LossTableManager::Instance()->FactorForAngleLimit()*CLHEP::hbarc/CLHEP::fermi;
  G4double p0 = 1.*hbarc/fermi;
//G4double cteta1 = 1.0 - p0*p0/2.0 * pow(A,2./3.)/(plab*plab);
  G4double cteta1 = 1.0 - p0*p0/2.0 * G4Pow::GetInstance()->Z23(A)/(plab*plab);
  if(cteta1 < -1.) cteta1 = -1.0;
  return cteta1;
}







