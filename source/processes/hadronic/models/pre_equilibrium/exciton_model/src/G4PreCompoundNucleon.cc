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
// by V. Lara
//
//J. M. Quesada (Apr. 2008) . Explicit inclusion of Coulomb barrier has been removed (NOW implicitely included through cross sections). Fixed bugs . New cross sections and combinatorial factor.

#include "G4PreCompoundNucleon.hh"
#include "G4PreCompoundParameters.hh"


G4double G4PreCompoundNucleon::
ProbabilityDistributionFunction(const G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) return 0.0;

//  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();

  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  


  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
 
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();


//JMQ (Jan. 08)fixed bugs
//  G4double A0 = (P*P+H*H+P-H)/4.0 - H/2.0;
  G4double A0 = ((P*P+H*H+P-H)/4.0 - H/2.0)/g0;
//  G4double A1 = A0 - P/2.0;
  G4double A1 = (A0 - P/2.0)/g1;
  
  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  G4double E1 = std::max(0.0,U - eKin - GetBindingEnergy() - A1);
  if (E1 == 0.0) return 0.0;

//JMQ commented for test30_04-06-08 OPT=3
//JMQ 02-06-08 Coulomb cutoff in xs is included here
//  if (eKin<theCoulombBarrier) return 0.0;

//JMQ (Dec.07) new cross sections added & combinatorial factor
 G4double Probability =std::pow(10.,-25.)* 1.0/pi2*2.0/(hbarc*hbarc*hbarc) * GetReducedMass()
* GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())  * eKin*CrossSection(eKin) * P*(N-1.0) * std::pow(g1*E1/(g0*E0),N-2.0)/E0 * g1/(g0*g0);

//JMQ & MAC (Jan. 08)  WARNING!!!!!!!!!!!!!
if (GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())<0.0 || CrossSection(eKin) <0) {
 G4cout<<"Rj="<<GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())
<<"   XS="<<CrossSection(eKin)<<G4endl;
G4cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<G4endl;}
/////



  return Probability;
}

////////////////////////////////////////////////////////////////////////////////////

//J. M. Quesada (Dec 2007-Apr 2008): New inverse reaction cross sections 
//OPT=1 Chatterjee's paramaterization for all ejectiles
//OPT=2     "                "         "  n,d,t,he3,alphas & Wellisch's parateterization for protons
//OPT=3 Kalbach's parameterization for all ejectiles
//OPT=4     "               "          "  n,d,t,he3,alphas & Wellisch's parateterization for protons
//

G4double G4PreCompoundNucleon::CrossSection(const  G4double K)
{
      G4double ResidualA=GetRestA();
      G4double ResidualZ=GetRestZ(); 
      G4int theA=static_cast<G4int>(GetA());
      G4int theZ=static_cast<G4int>(GetZ());
      G4double fragmentA=GetA()+GetRestA();

// OPT=2 //Default: Chatterjee's + Wellish's parameterizations
//JMQ 03-06-08 change to OPT=3 for n-Bi @ 63 MeV test
//JMQ 04-06-08 change to OPT=1 for test
      G4int OPT=1;

// Loop on XS options starts:
if ( OPT==1 ||OPT==2) { 
 
G4double Kc=K;

if (theA==1 && theZ==0) {
// Pramana (Bechetti & Greenles) for neutrons is chosen 

// JMQ  xsec is set constat above limit of validity
 if (K>50) Kc=50;

   G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs;

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
  else if (theA==1 && theZ==1 && OPT==1) {
 // Pramana ( Bechetti & Greenles) for protons is chosen

// JMQ  xsec is set constat above limit of validity
 if (K>50) Kc=50;

G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,xs;
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

G4double rnpro,rnneu,eekin,ekin,a,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th;
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
        if (xine_th < 0.) xine_th=0.0;        
        
if (xine_th <= 0.0 && eekin > 10.){
       std::ostringstream errOs;
   G4cout<<"WARNING: something funny happens with the Wellisch cross sections "<<G4endl; 
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}
        return xine_th;
}
      else {

        std::ostringstream errOs;
         G4cout<<"  Bad XS option in Pramana/Wellisch parameterizations"<<G4endl;
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());}
}

else if (OPT ==3 || OPT==4) {
//PRECO inverse cross sections are chosen

G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2;
G4double p, p0, p1, p2;

G4double flow,spill,xout,xpout,ares,athrd,ec,ecsq,xnulam,etest,ra,a,w,c,signor,signor2,sig; 
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
// PRECO xs for neutrons is choosen

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

    else if (theA==1 && theZ ==1 && OPT==3) {  
// PRECO xs for protons is choosen   

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
    
      ra = 0.;

                if (ares <= 60.)  signor = 0.92;
                else if (ares < 100.) signor = 0.8 + ares*0.002;


           xpout = GetZ();
  G4double  rz = ResidualZ;

           ec = 1.44 * xpout * rz / (1.5*athrd+ra);
           ecsq = ec * ec;
           p = p0 + p1/ec + p2/ecsq;
           landa = landa0*ares + landa1;
           a = std::pow(ares,mu1);
           mu = mu0 * a;
           nu = a* (nu0+nu1*ec+nu2*ecsq);
 
           c =std::min(3.15,ec*0.5);
           w = 0.7 * c / 3.15; 

           xnulam = nu / landa;
           if (xnulam > spill) xnulam=0.;
           if (xnulam < flow) goto loop1;

           etest =std::sqrt(xnulam) + 7.;

}  // end of  protons

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


   else {
        std::ostringstream errOs;
         G4cout<<"Bad ejectile in Kalbach parameterizations"<<G4endl;
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
        }
loop1:
      a = -2.*p*ec + landa - nu/ecsq;
      b = p*ecsq + mu + 2.*nu/ec;
      ecut = 0.;
      cut = a*a - 4.*p*b;
      if (cut > 0.) ecut = std::sqrt(cut);
      ecut = (ecut-a) / (p+p);
      ecut2 = ecut;
      if (cut < 0.) ecut2 = ecut - 2.;
      elab = K * fragmentA / ares;
      sig = 0.;
      if (elab <= ec) { //start for E<Ec 
           if (elab > ecut2)  sig = (p*elab*elab+a*elab+b) * signor;
             if(theA==1 && theZ==1) { //start for para protons
                     signor2 = (ec-elab-c) / w;
                     signor2 = 1. + std::exp(signor2);
                     sig = sig / signor2;
                                             }//end for protons      
            }              //end for E<Ec
      else {           //start for  E>Ec

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











