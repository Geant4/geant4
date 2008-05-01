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

#include "G4PreCompoundIon.hh"
#include "G4PreCompoundParameters.hh"


G4double G4PreCompoundIon::
ProbabilityDistributionFunction(const G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) return 0.0;
  
  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();

  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;

  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
 
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();



  G4double A0 = ((P*P+H*H+P-H)/4.0 - H/2.0)/g0; //JMQ fix

  G4double A1 = std::max(0.0,(A0*g0 + GetA()*(GetA()-2.0*P-1.0)/4.0)/g1); //JMQ fix

  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;

  G4double E1 = (std::max(0.0,GetMaximalKineticEnergy() - eKin - A1)); //JMQ re-fix 

  G4double Ej = std::max(0.0,eKin + GetBindingEnergy() ); //JMQ re-fix 


//JMQ combinatorial factors Rj and new cross sections have been added 
 G4double pA = std::pow(10.,-25.)*(3.0/4.0) * std::sqrt(std::max(0.0, 2.0/(GetReducedMass()*
(eKin+GetBindingEnergy()))))/(pi * r0 * r0 * std::pow(GetRestA(),2.0/3.0) )* eKin*CrossSection(eKin) /(r0*std::pow(GetRestA(),1.0/3.0)) * CoalescenceFactor(aFragment.GetA()) * FactorialFactor(N,P)* GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged())  ;

  G4double pB = std::pow((g1*E1)/(g0*E0),N-GetA()-1.0)*(g1/g0);

  G4double pC = std::pow((g1*Ej)/(g0*E0),GetA()-1.0)*(g1/g0)/E0; //JMQ fix

  G4double Probability = pA * pB * pC;

  return Probability;
}
////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-Apr 2008): New inverse reaction cross sections 
//OPT=1 Chatterjee's paramaterization for all ejectiles
//OPT=2     "                "         "  n,d,t,he3,alphas & Wellisch's parateterization for protons
//OPT=3 Kalbach's parameterization for all ejectiles

G4double G4PreCompoundIon::CrossSection(const G4double K)
{
      G4double ResidualA=GetRestA();
      G4double ResidualZ=GetRestZ(); 
      G4int theA=static_cast<G4int>(GetA());
      G4int theZ=static_cast<G4int>(GetZ());
      G4double fragmentA=GetA()+GetRestA();

// Default: Chatterjee's parameterization
      G4int OPT=3;

// Loop on XS options starts:
if ( OPT==1 ||OPT==2) {

G4double Kc=K;

// JMQ xsec is set constat above limit of validity
 if (K>50) Kc=50;


 G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,p, p0, p1, p2, Ec,delta,q,r,ji,xs;

 if (theA==2 && theZ==1) {
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


   if (xs <0.0) {xs=0.0;}
   return xs;
}
else if (OPT==3) {
//PRECO inverse cross sections are chosen

G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2,p, p0, p1, p2, signor,sig;
G4double flow,spill,xout,xpout,ares,athrd,ec,ecsq,xnulam,etest,ra,rz,a,w,c; 
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
 
 if (theA==2 && theZ==1) { // start deuterons
//     ** d from o.m. of perey and perey
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
           }   //end deuterons
 else if (theA==3 && theZ==1) { //start tritons

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
           }  //end tritons
else if (theA==3 && theZ==2) {//start he3



//c     ** 3he from o.m. of gibson et al
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
           }//end he3        
else if (theA==4 && theZ==2) {  //start  alphas
// c     ** alpha from huizenga and igo
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
           }  //end alphas

else  {

        std::ostringstream errOs;
         G4cout<<"  WARNING!! "<<G4endl;
        errOs << "theA=" << theA <<G4endl;
        errOs << "theZ=" << theZ <<G4endl;
        errOs << "this ejectile does not exist"  <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
      }
           ra=1.20;
           if (theA==2 || theA==3) ra=0.8; 
       
           xpout = theZ;
           rz = ResidualZ;
           ec = 1.44 * xpout * rz / (1.5*athrd+ra);
           ecsq = ec * ec;
           p = p0 + p1/ec + p2/ecsq;
           landa = landa0*ares + landa1;
           a = std::pow(ares,mu1);
           mu = mu0 * a;
           nu = a* (nu0+nu1*ec+nu2*ecsq);  
           xnulam = nu / landa;
           if (xnulam > spill) xnulam=0.;
           if (xnulam < flow) goto loop1;

           etest = 1.2 *std::sqrt(xnulam);

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
            }           //end for E<Ec
      else {           //start for E>Ec
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







