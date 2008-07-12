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
// $Id: G4EvaporationProbability.cc,v 1.15 2008-07-12 13:33:41 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
//J.M. Quesada (June 2008). Rebuilt class. Mayor changes: new inverse cross sections and numerical integration

#include "G4EvaporationProbability.hh"
#include "G4PairingCorrection.hh"
#include <iostream>

using namespace std;

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
//    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
//    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
//    G4double U = fragment.GetExcitationEnergy();
//    G4double CoulombBarrier = theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U);

    G4double MaximalKineticEnergy = anEnergy;
    G4double theEmissionProbability;

  if (MaximalKineticEnergy <= 0.0) 
    {
      theEmissionProbability = 0.0;
      return 0.0;
  }  


// JMQ June 08 LowerLimit is set again to 0 (Coulomb cutoff included in xs)
G4double LowerLimit=0.;
  

// JMQ: The MaximalKinetic energy just is over the barrier. Asimptotic value is needed.
//  G4double TrueMaximalKineticEnergy= MaximalKineticEnergy+CoulombBarrier;
//  already accounted for in G4EvaporationChannel (11/12/07)     

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
//JMQ (Dec. 2008) new method

G4double G4EvaporationProbability::ProbabilityDistributionFunction(const G4double K, const G4Fragment & fragment)
{ 

    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
  
    G4double U = fragment.GetExcitationEnergy();
   

   G4double delta01 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ));

   G4double delta00 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));


    G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));


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

//JMQ (June 08). Commented: Coulomb cutoff included, when needed,  in  xs  
//if(K<CoulombBarrier) {return 0.;}


//JMQ: 04-02-08 without 1/hbar_Panck remains as a width

Prob=std::pow(10.,-25.)*Gamma*NuclearMass/((pi*hbar_Planck)*(pi*hbar_Planck)*std::exp(2*std::sqrt(a*E0)))*K*CrossSection(K,fragment)*std::exp(2*std::sqrt(ap*E1));

return Prob;

}

////////////////////////////////////////////////////////////////////////////////////

//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=1 Chatterjee's paramaterization for all ejectiles+Coulomb cutoff
//OPT=2     "                "         "  n,d,t,he3,alphas & Wellisch's parateterization for protons + Coulomb cutoff
//OPT=3 Kalbach's parameterization for all ejectiles (improved for protons). NO Coulomb cutoff explicitely needed
//OPT=4     "               "          "  n,d,t,he3,alphas & Wellisch's parateterization for protons (Coulomb cutoff just for Wellisch's xs for protons)
// 
G4double G4EvaporationProbability::CrossSection(const G4double K,const  G4Fragment & fragment )
{

    G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
    G4double      athrd=std::pow(ResidualA,0.33333);

    G4double U = fragment.GetExcitationEnergy();
   
  G4double CoulombBarrier = theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U);

// Default OPT=2
      G4int OPT=2;


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

 
      landa = landa0/athrd + landa1;
      mu = mu0*athrd + mu1*athrd*athrd;
      nu = nu0*athrd*athrd*athrd*athrd + nu1*athrd*athrd + nu2 ;
       xs=landa*Kc + mu + nu/Kc;

if (xs < 0.0){
       std::ostringstream errOs;
   G4cout<<"WARNING:  NEGATIVE OPT=1 neutron cross section "<<G4endl;     
        errOs << "RESIDUAL: Ar=" << ResidualA << " Zr=" << ResidualZ <<G4endl;
        errOs <<"  xsec("<<Kc<<" MeV) ="<<xs <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
             }
         return xs;
           } 
  else if( theA==1 && theZ==1 && OPT==1) {

 // fitting to Bechetti & Greenles xs for protons is chosen

// JMQ 25/04/08 xsec is set constat above limit of validity
 if (K>50.) Kc=50.;
G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs(0.);
G4double p, p0, p1, p2,Ec,delta,q,r,ji;
//G4double Eo,epsilon1,epsilon2,discri;

        //JMQ (June 08) Coulomb cutoff
         if(K<=CoulombBarrier) return xs=0.0;   

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
      Ec = 1.44*theZ*ResidualZ/(1.5*athrd+delta);
      p = p0 + p1/Ec + p2/(Ec*Ec);
      landa = landa0*ResidualA + landa1;
      mu = mu0*std::pow(ResidualA,mu1);
      nu = std::pow(ResidualA,mu1)*(nu0 + nu1*Ec + nu2*Ec*Ec);
      q = landa - nu/(Ec*Ec) - 2*p*Ec;
      r = mu + 2*nu/Ec + p*Ec*Ec;



   ji=std::max(Kc,Ec);
   if(Kc < Ec) { xs = p*Kc*Kc + q*Kc + r;}
   else {xs = p*(Kc - ji)*(Kc - ji)+ landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}

  //JMQ 13-06-08 bug fixed unphysical of xs values removed
//JMQ 16-06-08 problems when Eo>Ec...commented and trivially solved 
//   discri=q*q-4.*p*r;
//   if(discri>=0) {
//   epsilon1=(-q+std::sqrt(discri))/2./p;
//   epsilon2=(-q-std::sqrt(discri))/2./p;
//   if(p>0.) Eo=std::max(epsilon1,epsilon2);
//   else    Eo=std::min(epsilon1,epsilon2);
//   if (Kc<Eo) return xs=0.;
//                 }
if (xs <0.0) {xs=0.0;}
//
//if (xs < 0.0 ){
//       std::ostringstream errOs;
//   G4cout<<"WARNING:  NEGATIVE OPT=1 proton cross section "<<G4endl; 
//        errOs << "RESIDUAL: Ar=" << ResidualA << " Zr=" << ResidualZ <<G4endl;
//        errOs <<"  xsec("<<Kc<<" MeV) ="<<xs <<G4endl;       
//        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
//              }
   return xs; 
}
   else if (theA==1 && theZ==1 && OPT==2) {
//Wellisch's  parameterization for protons is chosen
G4double rnpro,rnneu,eekin,ekin,a,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th(0.),athrdT;
 
        //JMQ June 08 Coulomb cutoff
         if(K<=CoulombBarrier) return xine_th=0.0;          

    eekin=K;
        rnpro=ResidualZ;
        rnneu=ResidualA-ResidualZ;
	a=rnneu+rnpro;
        athrdT=std::pow(a,0.333333);
	ekin=eekin/1000;
        r0=1.36*std::pow(static_cast<G4double>(10),static_cast<G4double>(-15));
	fac=pi*r0*r0;
	b0=2.247-0.915*(1.-1./athrdT);
	fac1=b0*(1.-1./athrdT);
	fac2=1.;
	if(rnneu > 1.5) fac2=std::log(rnneu);
	xine_th= std::pow(static_cast<G4double>(10),static_cast<G4double>(31))*fac*fac2*(1.+athrdT-fac1);
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
//        if (xine_th < 0.) {xine_th=0.0;} 
       if (xine_th < 0.0){
       std::ostringstream errOs;
       G4cout<<"WARNING: negative Wellisch cross section "<<G4endl; 
        errOs << "RESIDUAL: A=" << ResidualA << " Z=" << ResidualZ <<G4endl;
       errOs <<"  xsec("<<ekin<<" MeV) ="<<xine_th <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}  
     
        return xine_th;
} 
else {
G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs(0.);
G4double p, p0, p1, p2,Ec,delta,q,r,ji;
//G4double Eo,epsilon1,epsilon2,discri;
 
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

        //JMQ (June 08) Coulomb cutoff
         if(K<=CoulombBarrier) return xs=0.0;  
//
      Ec = 1.44*theZ*ResidualZ/(1.5*std::pow(ResidualA,0.333333)+delta);
      p = p0 + p1/Ec + p2/(Ec*Ec);
      landa = landa0*ResidualA + landa1;
      mu = mu0*std::pow(ResidualA,mu1);
      nu = std::pow(ResidualA,mu1)*(nu0 + nu1*Ec + nu2*(Ec*Ec));
      q = landa - nu/(Ec*Ec) - 2*p*Ec;
      r = mu + 2*nu/Ec + p*(Ec*Ec);

   ji=std::max(Kc,Ec);
   if(Kc < Ec) {xs = p*(Kc*Kc) + q*Kc + r;}
  else {xs = p*(Kc - ji)*(Kc - ji) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}


//JMQ 13-06-08 bug fixed unphysical negative xs values set to 0
//JMQ 16-06-08 problems when Eo>Ec...commented and trivially solved 
//   discri=q*q-4.*p*r;
//   if(discri>=0) {
//   epsilon1=(-q+std::sqrt(discri))/2./p;
//   epsilon2=(-q-std::sqrt(discri))/2./p;
//   if(p>0.) Eo=std::max(epsilon1,epsilon2);
//   else    Eo=std::min(epsilon1,epsilon2);
//   if (Kc<Eo) return xs=0.;                
//                 }
   if (xs <0.0) {xs=0.0;}
//
//if (xs < 0.0 ){
//       std::ostringstream errOs;
//   G4cout<<"WARNINGthe:  NEGATIVE OPT=1 || OPT=2 cluster cross section "<<G4endl; 
//        errOs << "EJECTILE: A=" << theA << " Z=" << theZ <<G4endl;
//        errOs << "RESIDUAL: Ar=" << ResidualA << " Zr=" << ResidualZ <<G4endl;
//        errOs <<"  xsec("<<Kc<<" MeV) ="<<xs <<G4endl;       
//        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
//                        }
   return xs;}
}
else if (OPT ==3 || OPT==4) {

G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2;
G4double p, p0, p1, p2;
G4double flow,spill,ec,ecsq,xnulam,etest,ra,a,w,c,signor,signor2,sig(0.); 
G4double b,ecut,cut,ecut2,geom,elab;
// JMQ (June 08) parameters for PRECO proton cross section refinement 
G4double afit,bfit,a2,b2;


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
afit=0.;
bfit=0.;
a2=0.;
b2=0.;
      flow = std::pow(10.,-18);
      spill= std::pow(10.,+18);
 
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

           if (ResidualA < 40.) signor=0.7+ResidualA*0.0075;
          if (ResidualA > 210.) signor = 1. + (ResidualA-210.)/250.;
          landa = landa0/athrd + landa1;
           mu = mu0*athrd + mu1*athrd*athrd;
           nu = nu0*athrd*ResidualA + nu1*athrd*athrd + nu2;
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

//JMQ 12/06/08  refinement of PRECO proton cross sections
//            afit=0.00464447;
//            bfit=2.94443;
//            a2= -0.00211134;
//            b2= 0.0548949;    
//JMQ 15/06/08  
            afit=-0.0785656;
            bfit=5.10789;
            a2= -0.00089076;
            b2= 0.0231597;  
         } //end protons with PRECO's param.

 else if (theA==1 && theZ ==1 && OPT==4) { //Wellisch's  parameterization for protons is chosen
 G4double rnpro,rnneu,eekin,ekin,a,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th(0.),athrdT;

        //JMQ (June 08) Coulomb cutoff
        if(K<=CoulombBarrier) return xine_th=0.0;          

        eekin=K;
        rnpro=ResidualZ;
        rnneu=ResidualA-ResidualZ;
	a=rnneu+rnpro;
        athrdT=std::pow(a,0.33333);
	ekin=eekin/1000;
        r0=1.36*std::pow(static_cast<G4double>(10),static_cast<G4double>(-15));
	fac=pi*std::pow(r0,2);
	b0=2.247-0.915*(1.-1./athrdT);
	fac1=b0*(1.-1./athrdT);
	fac2=1.;
	if(rnneu > 1.5) fac2=std::log(rnneu);
	xine_th= std::pow(static_cast<G4double>(10),static_cast<G4double>(31))*fac*fac2*(1.+athrdT-fac1);
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
//        if (xine_th < 0.) {xine_th=0.0;}    
if (xine_th < 0.0){
       std::ostringstream errOs;
   G4cout<<"WARNING: negative Wellisch cross section "<<G4endl; 
        errOs << "RESIDUAL: A=" << ResidualA << " Z=" << ResidualZ <<G4endl;
       errOs <<"  xsec("<<ekin<<" MeV) ="<<xine_th <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}

        return xine_th;     
}
 
 else if (theA==2 && theZ==1) { //start deuterons

        //JMQ (June 08) Coulomb cutoff
 //        if(K<=CoulombBarrier) return sig=0.0;  

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

        //JMQ (June 08) Coulomb cutoff
 //        if(K<=CoulombBarrier) return sig=0.0;  

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

        //JMQ (June 08) Coulomb cutoff
//         if(K<=CoulombBarrier) return sig=0.0;  

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

        //JMQ (June 08) Coulomb cutoff
//         if(K<=CoulombBarrier) return sig=0.0;  

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
      errOs << "WARNING: OPT=3 || OPT=4: this ejectile does not exist"  <<G4endl;
     throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
 
//still inside charged particles loop                   
                 ra=1.20;
//just for protons
           if (theA==1 && theZ==1) {
                 ra = 0.;
                if (ResidualA <= 60.)  signor = 0.92;
                else if (ResidualA < 100.) signor = 0.8 + ResidualA*0.002;
          }

           ec = 1.44 * theZ * ResidualZ / (1.5*athrd+ra);
           ecsq = ec * ec;
           p = p0 + p1/ec + p2/ecsq;
           landa = landa0*ResidualA + landa1;
           a = std::pow(ResidualA,mu1);
          mu = mu0 * a;
           nu = a* (nu0+nu1*ec+nu2*ecsq);
           if (theA==2 || theA==3) ra=0.8;
             if (theA==1 && theZ==1) {
                c =std::min(3.15,ec*0.5);
                w = 0.7 * c / 3.15; }
           xnulam = nu / landa;
           if (xnulam > spill) xnulam=0.;
           if (xnulam < flow) goto loop1;
           if(theA==1 && theZ==1)  etest =std::sqrt(xnulam) + 7.;
 
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
     elab = K * fragment.GetA() / ResidualA;
      sig = 0.;

      if (elab <= ec) {  
           if (elab > ecut2)  sig = (p*elab*elab+a*elab+b) * signor;
             if(theA==1 && theZ==1) { 
                     signor2 = (ec-elab-c) / w;
                     signor2 = 1. + std::exp(signor2);
                     sig = sig / signor2; 
// signor2 is empirical global corr'n at low elab for protons in PRECO, not enough for p+27Al
//JMQ 12/06/08  refinement for proton cross section
   if (ResidualZ<=26) {sig = sig*std::exp(-(a2*ResidualZ + b2)*(elab-(afit*ResidualZ+bfit)*ec)*(elab-(afit*ResidualZ+bfit)*ec)); 
}
                                      } 
                       }              
      else {sig = (landa*elab+mu+nu/elab) * signor;
//JMQ 12/06/08  refinement for proton cross section
   if (theA==1 && theZ==1 && ResidualZ<=26 && elab <=(afit*ResidualZ+bfit)*ec) {sig = sig*std::exp(-(a2*ResidualZ + b2)*(elab-(afit*ResidualZ+bfit)*ec)*(elab-(afit*ResidualZ+bfit)*ec));}
//
           geom = 0.;
         if (xnulam < flow || elab < etest) return sig;
           geom = std::sqrt(theA*K);
           geom = 1.23*athrd + ra + 4.573/geom;
           geom = 31.416 * geom * geom;
           sig = std::max(geom,sig);
      }
return sig;} 
else 
{
std::ostringstream errOs;
      errOs << "BAD CROSS SECTION OPTION !!"  <<G4endl;
     throw G4HadronicException(__FILE__, __LINE__, errOs.str());
return 0.;
}
}
