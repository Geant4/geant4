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
// $Id: G4EvaporationProbability.cc,v 1.12 2008-06-05 16:46:29 quesada Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// JMQ & MAC 07/12/2007: New inverse cross sections
//
//J.M. Quesada (Dec 2007-Apr 2008). Rebuilt class. Mayor changes: new inverse cross sections and 
//numerical integration

#include <iostream>
using namespace std;

#include "G4EvaporationProbability.hh"
#include "G4PairingCorrection.hh"



G4EvaporationProbability::G4EvaporationProbability(const G4EvaporationProbability &) : G4VEmissionProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationProbability::copy_constructor meant to not be accessable");
}




const G4EvaporationProbability & G4EvaporationProbability::
operator=(const G4EvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4EvaporationProbability::operator==(const G4EvaporationProbability &) const
{
    return false;
}

G4bool G4EvaporationProbability::operator!=(const G4EvaporationProbability &) const
{
    return true;
}
  
G4double G4EvaporationProbability::EmissionProbability(const G4Fragment & fragment, const G4double anEnergy)
{
    G4double EmissionProbability = 0.0;
    G4double MaximalKineticEnergy = anEnergy;

    if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
	EmissionProbability = CalculateProbability(MaximalKineticEnergy,fragment);

    }
    return EmissionProbability;
}


//////////////////////////////////////////////////////////////////////////////
//JMQ: rebuilt method

G4double G4EvaporationProbability::CalculateProbability(const G4double anEnergy,const G4Fragment & fragment)
{

    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);

    G4double U = fragment.GetExcitationEnergy();
    G4double CoulombBarrier = theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U);

    G4double MaximalKineticEnergy = anEnergy;
    G4double theEmissionProbability;

  if (MaximalKineticEnergy <= 0.0) 
    {
      theEmissionProbability = 0.0;
      return 0.0;
  }  
// JMQ:  
// Now: integration is numericaly performed & cross section is set to 0. when negative =>
// JMQ May.08 LowerLimit is set again to Coulomb Barrier
  G4double LowerLimit = CoulombBarrier;

  

// JMQ: The MaximalKinetic energy just is over the barrier. Asimptotic value is needed.
//  G4double TrueMaximalKineticEnergy= MaximalKineticEnergy+CoulombBarrier;
//  already accounted  in G4EvaporationChannel (11/12/07)     

  G4double UpperLimit = MaximalKineticEnergy;


  G4double Width = IntegrateEmissionProbability(LowerLimit,UpperLimit,fragment);
  
  return Width;
  
}

/////////////////////////////////////////////////////////////////////

G4double G4EvaporationProbability::
IntegrateEmissionProbability(const G4double & Low, const G4double & Up, const G4Fragment & fragment)
{

  static const G4int N = 10;
  // 10-Points Gauss-Legendre abcisas and weights
  static const G4double w[N] = {
    0.0666713443086881,
    0.149451349150581,
    0.219086362515982,
    0.269266719309996,
    0.295524224714753,
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.0666713443086881
  };
  static const G4double x[N] = {
    -0.973906528517172,
    -0.865063366688985,
    -0.679409568299024,
    -0.433395394129247,
    -0.148874338981631,
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172
  };

  G4double Total = 0.0;


  for (G4int i = 0; i < N; i++) 
    {

      G4double KineticE = ((Up-Low)*x[i]+(Up+Low))/2.0;

      Total += w[i]*ProbabilityDistributionFunction(KineticE,fragment);

    }
  Total *= (Up-Low)/2.0;


  return Total;
}


///////////////////////////////////////////////////////////////////////
//JMQ Dec. 2008 new method

G4double G4EvaporationProbability::ProbabilityDistributionFunction(const G4double K, const G4Fragment & fragment)
{ 

    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
  
    G4double U = fragment.GetExcitationEnergy();
   
 G4double CoulombBarrier = theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U);


   G4double delta01 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ));

   G4double delta00 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));

//
    G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));

//    G4double U = fragment.GetExcitationEnergy();

   G4double NuclearMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass(theZ,theA);

   G4double theSeparationEnergy= G4NucleiProperties::GetMassExcess(static_cast<G4int>(theA),static_cast<G4int>(theZ)) +
    G4NucleiProperties::GetMassExcess(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ)) -
     G4NucleiProperties::GetMassExcess(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));

   G4double ap = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ),U - theSeparationEnergy - delta0);

   G4double a = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()),U - delta0);

G4double Prob;



G4double E1=U-theSeparationEnergy-delta01-K;

if (E1<0.) return 0.;

G4double E0=U-delta00;

//JMQ commented for test30_04-06-08 OPT=3
//JMQ: 04-05-08 cross sections are set to 0 below Coulomb Barrier

//if(K<CoulombBarrier) {return 0.;}



//JMQ: 04-02-08 without 1/hbar_Panck remains as a width

Prob=std::pow(10.,-25.)*Gamma*NuclearMass/(std::pow(pi*hbar_Planck,2.)*std::exp(2*std::sqrt(a*E0)))*K*CrossSection(K,fragment)*std::exp(2*std::sqrt(ap*E1));

return Prob;

}

////////////////////////////////////////////////////////////////////////////////////

//J. M. Quesada (Dec 2007-Apr 2008): New inverse reaction cross sections 
//OPT=1 Chatterjee's paramaterization for all ejectiles
//OPT=2     "                "         "  n,d,t,he3,alphas & Wellisch's parateterization for protons
//OPT=3 Kalbach's parameterization for all ejectiles
//OPT=4     "               "          "  n,d,t,he3,alphas & Wellisch's parateterization for protons
// 
G4double G4EvaporationProbability::CrossSection(const G4double K,const  G4Fragment & fragment )
{

    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
    G4int theA=static_cast<G4int>(GetA());
    G4int theZ=static_cast<G4int>(GetZ());


// OPT=2 //Default: Chatterjee's + Wellish's parameterizations
//JMQ 03-06-08 change to OPT=3 for n-Bi @ 63 MeV test
//JMQ 04-06-08 change to OPT=1 for test
      G4int OPT=1;

// Loop on XS options starts:
if ( OPT==1 ||OPT==2) { 

G4double Kc=K;

  if (theA==1  && theZ==0) {
// fitting to Bechetti & Greenles xs  for neutrons is chosen 

// JMQ 25/04/08 xsec is set constat above limit of validity
 if (K>50) Kc=50;

   G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs(0.);

       landa0 = 18.57;
       landa1 = -22.93;
       mu0 = 381.7;
       mu1 = 24.31;
       nu0 = 0.172;
       nu1 = -15.39;
       nu2 = 804.8; 
      landa = landa0*std::pow(ResidualA,-0.33333) + landa1;
      mu = mu0*std::pow(ResidualA,0.333333) + mu1*std::pow(ResidualA,0.666666);
      nu = nu0*std::pow(ResidualA,4./3.) + nu1*std::pow(ResidualA,0.666666) + nu2 ;
       xs=landa*Kc + mu + nu/Kc;

if (xs <= 0.0 && Kc > 10){
       std::ostringstream errOs;
   G4cout<<"WARNING: something funny happens with the Pramana neutron cross section "<<G4endl; 
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}
         return xs;
           } 
  else if( theA==1 && theZ==1 && OPT==1) {

 // fitting to Bechetti & Greenles xs for protons is chosen

// JMQ 25/04/08 xsec is set constat above limit of validity
 if (K>50.) Kc=50.;
G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs(0.);
G4double p, p0, p1, p2,Ec,delta,q,r,ji;
       p0 = 15.72;
      p1 = 9.65;
      p2 = -449.0;
      landa0 = 0.00437;
      landa1 = -16.58;
      mu0 = 244.7;
      mu1 = 0.503;
      nu0 = 273.1;
      nu1 = -182.4;
      nu2 = -1.872;  
      delta=0.;  
      Ec = 1.44*theZ*ResidualZ/(1.5*std::pow(ResidualA,0.333333)+delta);
      p = p0 + p1/Ec + p2/std::pow(Ec,2.);
      landa = landa0*ResidualA + landa1;
      mu = mu0*std::pow(ResidualA,mu1);
      nu = std::pow(ResidualA,mu1)*(nu0 + nu1*Ec + nu2*std::pow(Ec,2.));
      q = landa - nu/std::pow(Ec,2.) - 2*p*Ec;
      r = mu + 2*nu/Ec + p*std::pow(Ec,2.);
   ji=std::max(Kc,Ec);
   if(Kc < Ec) { xs = p*std::pow(Kc,2.) + q*Kc + r;}
   else {xs = p*std::pow((Kc - ji),2.) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}

  //JMQ 05-06-08 bug fixed unphysical of xs values removed
   G4double Eo,epsilon1,epsilon2;
   epsilon1=(-q+std::sqrt(q*q-4.*p*r))/2./p;
   epsilon2=(-q-std::sqrt(q*q-4.*p*r))/2./p;
   if(p>0.) Eo=std::max(epsilon1,epsilon2);
   else    Eo=std::min(epsilon1,epsilon2);
   if (Kc<Eo) xs=0.;
 //

   if (xs <0.0) {xs=0.0;}
if (xs <= 0.0 && Kc > 2*Ec){
       std::ostringstream errOs;
   G4cout<<"WARNING: something funny happens with the Pramana proton cross section "<<G4endl; 
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}
   return xs; 
}
   else if (theA==1 && theZ==1 && OPT==2) {
//Wellisch's  parameterization for protons is chosen
G4double rnpro,rnneu,eekin,ekin,a,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th(0.);
    eekin=K;
        rnpro=ResidualZ;
        rnneu=ResidualA-ResidualZ;
	a=rnneu+rnpro;
	ekin=eekin/1000;
        r0=1.36*std::pow(static_cast<G4double>(10),static_cast<G4double>(-15));
	fac=pi*std::pow(r0,2);
	b0=2.247-0.915*(1.-std::pow(a,-0.33333));
	fac1=b0*(1.-std::pow(a,-0.333333));
	fac2=1.;
	if(rnneu > 1.5) fac2=std::log(rnneu);
	xine_th= std::pow(static_cast<G4double>(10),static_cast<G4double>(31))*fac*fac2*(1.+std::pow(a,0.3333)-fac1);
       xine_th=(1.-0.15*std::exp(-ekin))*xine_th/(1.00-0.0007*a);	
       ff1=0.70-0.0020*a ;
       ff2=1.00+1/a;
       ff3=0.8+18/a-0.002*a;
       fac=1.-(1./(1.+std::exp(-8.*ff1*(std::log10(ekin)+1.37*ff2))));
       xine_th=xine_th*(1.+ff3*fac);
       	ff1=1.-1/a-0.001*a;
	ff2=1.17-2.7/a-0.0014*a;
	fac=-8.*ff1*(std::log10(ekin)+2.0*ff2);
	fac=1./(1.+std::exp(fac));
	xine_th=xine_th*fac;
        if (xine_th < 0.) {xine_th=0.0;}        
        return xine_th;
} 
else {
G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs(0.);
G4double p, p0, p1, p2,Ec,delta,q,r,ji;
 
  if (theA==2 && theZ==1) {
// JMQ 25/04/08 parameterization limit
 if (K>50.) Kc=50.;
       p0 = -38.21;
      p1 = 922.6;
      p2 = -2804.;
      landa0 = -0.0323;
      landa1 = -5.48;
      mu0 = 336.1;
      mu1 = 0.48;
      nu0 = 524.3;
      nu1 = -371.8;
      nu2 = -5.924;  
      delta=1.2;            
           }     
 else if (theA==3 && theZ==1) {
// JMQ 25/04/08 parameterization limit
 if (K>50.) Kc=50.;
        p0 = -11.04;
      p1 = 619.1;
      p2 = -2147.;
      landa0 = -0.0426;
      landa1 = -10.33;
      mu0 = 601.9;
      mu1 = 0.37;
      nu0 = 583.0;
      nu1 = -546.2;
      nu2 = 1.718;  
      delta=1.2;               
           }    
 else if (theA==3 && theZ==2) {
// JMQ 25/04/08 parametarization limit
 if (K>50.) Kc=50.;
        p0 = -3.06;
      p1 = 278.5;
      p2 = -1389.;
      landa0 = -0.00535;
      landa1 = -11.16;
      mu0 = 555.5;
      mu1 = 0.40;
      nu0 = 687.4;
      nu1 = -476.3;
      nu2 = 0.509;    
      delta=1.2;             
           }
 else if (theA==4 && theZ==2) {
// JMQ 25/04/08 parameterization limit
 if (K>50.) Kc=50.;

         p0 = 10.95;
      p1 = -85.2;
      p2 = 1146.;
      landa0 = 0.0643;
      landa1 = -13.96;
      mu0 = 781.2;
      mu1 = 0.29;
      nu0 = -304.7;
      nu1 = -470.0;
      nu2 = -8.580;   
      delta=1.2;            
           }
else  {  
     std::ostringstream errOs;
      errOs << "theA=" << theA <<G4endl;
      errOs << "theZ=" << theZ <<G4endl;
      errOs << "this ejectile does not exist"  <<G4endl;
     throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
      Ec = 1.44*theZ*ResidualZ/(1.5*std::pow(ResidualA,0.333333)+delta);
      p = p0 + p1/Ec + p2/std::pow(Ec,2.);
      landa = landa0*ResidualA + landa1;
      mu = mu0*std::pow(ResidualA,mu1);
      nu = std::pow(ResidualA,mu1)*(nu0 + nu1*Ec + nu2*std::pow(Ec,2.));
      q = landa - nu/std::pow(Ec,2.) - 2*p*Ec;
      r = mu + 2*nu/Ec + p*std::pow(Ec,2.);
   ji=std::max(Kc,Ec);
   if(Kc < Ec) { xs = p*std::pow(Kc,2.) + q*Kc + r;}
   else {xs = p*std::pow((Kc - ji),2.) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}


  //JMQ 05-06-08 bug fixed unphysical of xs values removed
   G4double Eo,epsilon1,epsilon2;
   epsilon1=(-q+std::sqrt(q*q-4.*p*r))/2./p;
   epsilon2=(-q-std::sqrt(q*q-4.*p*r))/2./p;
   if(p>0.) Eo=std::max(epsilon1,epsilon2);
   else    Eo=std::min(epsilon1,epsilon2);
   if (Kc<Eo) xs=0.;
 //

   if (xs <0.0) {xs=0.0;}
   return xs;}
}
else if (OPT ==3 || OPT==4) {

G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2;
G4double p, p0, p1, p2;

G4double flow,spill,xout,xpout,ares,athrd,ec,ecsq,xnulam,etest,ra,a,w,c,signor,signor2,sig(0.); 
G4double b,ecut,cut,ecut2,geom,elab;
G4int jout,jpout,jnout;

//safety initialization
landa0=0;
landa1=0;
mu0=0.;
mu1=0.;
nu0=0.;
nu1=0.;
nu2=0.;
p0=0.;
p1=0.;
p2=0.;
etest=0.;
w=0.;
c=0.;
      flow = std::pow(10.,-18);
      spill= std::pow(10.,+18);

      jout=theA;
      jpout=theZ;
      jnout=jout-jpout;
      xout=GetA();
      ares=ResidualA;
      athrd = std::pow(ares,0.3333);

      signor = 1.;

   if (theA==1  && theZ==0) {
// PRECO neutron are choosen

      p0 = -312.;
      p1= 0.;
      p2 = 0.;
      landa0 = 12.10;
      landa1=  -11.27;
      mu0 = 234.1;
      mu1 = 38.26;
      nu0 = 1.55;
      nu1 = -106.1;
      nu2 = 1280.8; 

           if (ares < 40.) signor=0.7+ares*0.0075;
          if (ares > 210.) signor = 1. + (ares-210.)/250.;
          landa = landa0/athrd + landa1;
           mu = mu0*athrd + mu1*athrd*athrd;
           nu = nu0*athrd*ares + nu1*athrd*athrd + nu2;
           ec = 0.5;
           ecsq = 0.25;
          p = p0;
           xnulam = 1.;
           etest = 32.;
//          ** etest is the energy above which the rxn cross section is
//          ** compared with the geometrical limit and the max taken.
//          ** xnulam here is a dummy value to be used later.

          ra=0.;
     
} // end if of PRECO neutrons

    else {
     if (theA==1 && theZ ==1 && OPT==3) {  
// PRECO Protons are choosen


//     ** p from  becchetti and greenlees (but modified with sub-barrier
//     ** correction function and xp2 changed from -449)

     p0 = 15.72;
      p1 = 9.65;
      p2 = -300.;
     landa0 = 0.00437;
      landa1 = -16.58;
      mu0 = 244.7;
      mu1 = 0.503;
      nu0 = 273.1;
     nu1 = -182.4;
      nu2 = -1.872;

         } //end protons with PRECO's param.

 else if (theA==1 && theZ ==1 && OPT==4) { //Wellisch's  parameterization for protons is chosen
 G4double rnpro,rnneu,eekin,ekin,a,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th(0.);
    eekin=K;
        rnpro=ResidualZ;
        rnneu=ResidualA-ResidualZ;
	a=rnneu+rnpro;
	ekin=eekin/1000;
        r0=1.36*std::pow(static_cast<G4double>(10),static_cast<G4double>(-15));
	fac=pi*std::pow(r0,2);
	b0=2.247-0.915*(1.-std::pow(a,-0.33333));
	fac1=b0*(1.-std::pow(a,-0.333333));
	fac2=1.;
	if(rnneu > 1.5) fac2=std::log(rnneu);
	xine_th= std::pow(static_cast<G4double>(10),static_cast<G4double>(31))*fac*fac2*(1.+std::pow(a,0.3333)-fac1);
       xine_th=(1.-0.15*std::exp(-ekin))*xine_th/(1.00-0.0007*a);	
       ff1=0.70-0.0020*a ;
       ff2=1.00+1/a;
       ff3=0.8+18/a-0.002*a;
       fac=1.-(1./(1.+std::exp(-8.*ff1*(std::log10(ekin)+1.37*ff2))));
       xine_th=xine_th*(1.+ff3*fac);
       	ff1=1.-1/a-0.001*a;
	ff2=1.17-2.7/a-0.0014*a;
	fac=-8.*ff1*(std::log10(ekin)+2.0*ff2);
	fac=1./(1.+std::exp(fac));
	xine_th=xine_th*fac;
        if (xine_th < 0.) {xine_th=0.0;}        
        return xine_th;
}
 
 else if (theA==2 && theZ==1) { //start deuterons

     p0 = 0.798;
      p1 = 420.3;
     p2 = -1651.;
      landa0 = 0.00619;
      landa1 = -7.54;
      mu0 = 583.5;
      mu1 = 0.337;
      nu0 = 421.8;
      nu1 = -474.5;
     nu2 = -3.592;      
           }                           //end deuterons
 else if (theA==3 && theZ==1) {       //start  tritons

//     ** t from o.m. of hafele, flynn et al
      p0 = -21.45;
    p1 = 484.7;
     p2 = -1608.;
      landa0 = 0.0186;
      landa1 = -8.90;
      mu0 = 686.3;
      mu1 = 0.325;
      nu0 = 368.9;
      nu1 = -522.2;
      nu2 = -4.998;

           }                     //end tritons

else if (theA==3 && theZ==2) {   //start he3

//     ** 3he from o.m. of gibson et al
     p0 = -2.88;
      p1 = 205.6;
     p2 = -1487.;
     landa0 = 0.00459;
      landa1 = -8.93;
      mu0 = 611.2;
      mu1 = 0.35;
      nu0 = 473.8;
     nu1 = -468.2;
      nu2 = -2.225;        
           }               //end he3 he3
      
  
else if (theA==4 && theZ==2) { //start  alphas

//      ** alpha from huizenga and igo
        p0 = 10.95;
      p1 = -85.2;
      p2 = 1146.;
      landa0 = 0.0643;
      landa1 = -13.96;
      mu0 = 781.2;
      mu1 = 0.29;
      nu0 = -304.7;
      nu1 = -470.0;
      nu2 = -8.580;        
           }  //end  alphas
else  {  
     std::ostringstream errOs;
      errOs << "theA=" << theA <<G4endl;
      errOs << "theZ=" << theZ <<G4endl;
      errOs << "this ejectile does not exist"  <<G4endl;
     throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
 
//still inside charged particles loop                   
                 ra=1.20;
//just for protons
           if (theA==1 && theZ==1) {
                 ra = 0.;
                if (ares <= 60.)  signor = 0.92;
                else if (ares < 100.) signor = 0.8 + ares*0.002;
          }
           xpout = theZ;
 G4double  rz = ResidualZ;

           ec = 1.44 * xpout * rz / (1.5*athrd+ra);
           ecsq = ec * ec;
           p = p0 + p1/ec + p2/ecsq;
           landa = landa0*ares + landa1;
           a = std::pow(ares,mu1);
          mu = mu0 * a;
           nu = a* (nu0+nu1*ec+nu2*ecsq);
           if (jout==2 || jout==3) ra=0.8;
             if (jpout==1 && jnout==0) {
                c =std::min(3.15,ec*0.5);
                w = 0.7 * c / 3.15; }
           xnulam = nu / landa;
           if (xnulam > spill) xnulam=0.;
           if (xnulam < flow) goto loop1;
           if(jpout==1 && jnout==0)  etest =std::sqrt(xnulam) + 7.;
 
// and for the rest of charged
            else    etest = 1.2 *std::sqrt(xnulam);

//          ** For xnulam.gt.0, sig reaches a maximum at sqrt(xnulam).

};  // end of else of charged

loop1:
      a = -2.*p*ec + landa - nu/ecsq;
      b = p*ecsq + mu + 2.*nu/ec;
      ecut = 0.;
      cut = a*a - 4.*p*b;
      if (cut > 0.) ecut = std::sqrt(cut);
      ecut = (ecut-a) / (p+p);
      ecut2 = ecut;
      if (cut < 0.) ecut2 = ecut - 2.;
     elab = K * fragment.GetA() / ares;
      sig = 0.;

      if (elab <= ec) {  
           if (elab > ecut2)  sig = (p*elab*elab+a*elab+b) * signor;
             if(jpout==1 && jnout==0) { 
                     signor2 = (ec-elab-c) / w;
                     signor2 = 1. + std::exp(signor2);
                     sig = sig / signor2;                     }   
            }              
      else {           
           sig = (landa*elab+mu+nu/elab) * signor;
           geom = 0.;
//           if (xnulam < flow) goto loop2; 
//           if (elab < etest) goto loop2;
         if (xnulam < flow || elab < etest) return sig;
           geom = std::sqrt(xout*K);
           geom = 1.23*athrd + ra + 4.573/geom;
           geom = 31.416 * geom * geom;
           sig = std::max(geom,sig);
      }

//loop2:
 return sig;} 
else 
{
std::ostringstream errOs;
      errOs << "BAD CROSS SECTION OPTION !!"  <<G4endl;
     throw G4HadronicException(__FILE__, __LINE__, errOs.str());
return 0.;
}
}
