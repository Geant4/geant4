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
//J.M.Quesada (August 08). New source file
// 
// Modif (21 August 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// 
// JMQ (06 September 2008) Also external choice has been added for:
//                      - superimposed Coulomb barrier (if useSICB=true) 

#include "G4PreCompoundProton.hh"

  G4ReactionProduct * G4PreCompoundProton::GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct = 
      new G4ReactionProduct(G4Proton::ProtonDefinition());
    theReactionProduct->SetMomentum(GetMomentum().vect());
    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef PRECOMPOUND_TEST
    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
    return theReactionProduct;
  }


   G4double G4PreCompoundProton::GetRj(const G4int NumberParticles, const G4int NumberCharged)
  {
    G4double rj = 0.0;
    if(NumberParticles > 0) rj = static_cast<G4double>(NumberCharged)/static_cast<G4double>(NumberParticles);
    return rj;
  }



////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1 Chatterjee's paramaterization 
//OPT=2 Wellisch's parametarization
//OPT=3 Kalbach's parameterization 
// 
 G4double G4PreCompoundProton::CrossSection(const  G4double K)
{
  //G4cout<<" In G4PreCompoundProton OPTxs="<<OPTxs<<G4endl;
  //G4cout<<" In G4PreCompoundProton useSICB="<<useSICB<<G4endl;

  ResidualA=GetRestA();
  ResidualZ=GetRestZ(); 
  theA=GetA();
  theZ=GetZ();
  ResidualAthrd=std::pow(ResidualA,0.33333);
  FragmentA=GetA()+GetRestA();
  FragmentAthrd=std::pow(FragmentA,0.33333);



  if (OPTxs==0) return GetOpt0(K);
  else if( OPTxs==1) return GetOpt1(K);
  else if( OPTxs==2|| OPTxs==4) return GetOpt2(K);
  else if (OPTxs==3)  return GetOpt3(K);
  else{
    std::ostringstream errOs;
    errOs << "BAD PROTON CROSS SECTION OPTION !!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
}

// *********************** OPT=0 : Dostrovski's cross section  *****************************

G4double G4PreCompoundProton::GetOpt0(const  G4double K)
{
const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  //with the rest of the options
 return 1.e+25*pi*(r0*ResidualAthrd)*(r0*ResidualAthrd)*GetAlpha()*(1.+GetBeta()/K);
}
//
//------------
//
  G4double G4PreCompoundProton::GetAlpha()
  {
 G4double aZ = static_cast<G4double>(GetRestZ());
    G4double C = 0.0;
    if (aZ >= 70) 
      {
	C = 0.10;
      } 
    else 
      {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
      }
    return 1.0 + C;
  }
//
//-------------------
//  
  G4double G4PreCompoundProton::GetBeta() 
  {
  return -GetCoulombBarrier();
  }
//
  


//********************* OPT=1 : Chatterjee's cross section ************************ 
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundProton::GetOpt1(const  G4double K)
{
 G4double Kc=K; 

// JMQ  xsec is set constat above limit of validity
 if (K>50)  Kc=50;

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

 Ec = 1.44*theZ*ResidualZ/(1.5*ResidualAthrd+delta);
 p = p0 + p1/Ec + p2/(Ec*Ec);
 landa = landa0*ResidualA + landa1;
 mu = mu0*std::pow(ResidualA,mu1);
 nu = std::pow(ResidualA,mu1)*(nu0 + nu1*Ec + nu2*(Ec*Ec));
 q = landa - nu/(Ec*Ec) - 2*p*Ec;
 r = mu + 2*nu/Ec + p*(Ec*Ec);

 ji=std::max(Kc,Ec);
 if(Kc < Ec) { xs = p*Kc*Kc + q*Kc + r;}
 else {xs = p*(Kc - ji)*(Kc - ji) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}
 if (xs <0.0) {xs=0.0;}

 return xs; 

}



//************* OPT=2 : Welisch's proton reaction cross section ************************ 

G4double G4PreCompoundProton::GetOpt2(const  G4double K)
{

  G4double rnpro,rnneu,eekin,ekin,ff1,ff2,ff3,r0,fac,fac1,fac2,b0,xine_th(0);
 
//This is redundant when the Coulomb  barrier is overimposed to all cross sections 
//It should be kept when Coulomb barrier only imposed at OPTxs=2

  if(!useSICB && K<=theCoulombBarrier) return xine_th=0.0;

  eekin=K;
  rnpro=ResidualZ;
  rnneu=ResidualA-ResidualZ;
  ekin=eekin/1000;
  r0=1.36*1.e-15;
  fac=pi*r0*r0;
  b0=2.247-0.915*(1.-1./ResidualAthrd);
  fac1=b0*(1.-1./ResidualAthrd);
  fac2=1.;
  if(rnneu > 1.5) fac2=std::log(rnneu);
  xine_th= 1.e+31*fac*fac2*(1.+ResidualAthrd-fac1);
  xine_th=(1.-0.15*std::exp(-ekin))*xine_th/(1.00-0.0007*ResidualA);	
  ff1=0.70-0.0020*ResidualA ;
  ff2=1.00+1/ResidualA;
  ff3=0.8+18/ResidualA-0.002*ResidualA;
  fac=1.-(1./(1.+std::exp(-8.*ff1*(std::log10(ekin)+1.37*ff2))));
  xine_th=xine_th*(1.+ff3*fac);
  ff1=1.-1/ResidualA-0.001*ResidualA;
  ff2=1.17-2.7/ResidualA-0.0014*ResidualA;
  fac=-8.*ff1*(std::log10(ekin)+2.0*ff2);
  fac=1./(1.+std::exp(fac));
  xine_th=xine_th*fac;            
  if (xine_th < 0.0){
    std::ostringstream errOs;
    G4cout<<"WARNING:  negative Wellisch cross section "<<G4endl; 
    errOs << "RESIDUAL: A=" << ResidualA << " Z=" << ResidualZ <<G4endl;
    errOs <<"  xsec("<<ekin<<" MeV) ="<<xine_th <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
  return xine_th;
            
}


// *********** OPT=3 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundProton::GetOpt3(const  G4double K)
{
  //     ** p from  becchetti and greenlees (but modified with sub-barrier
  //     ** correction function and xp2 changed from -449)
  // JMQ (june 2008) : refinement of proton cross section for light systems
  //
  G4double landa, landa0, landa1, mu, mu0, mu1,nu, nu0, nu1, nu2;
  G4double p, p0, p1, p2;
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
  
  // parameters for  proton cross section refinement 
  G4double afit,bfit,a2,b2;
  afit=-0.0785656;
  bfit=5.10789;
  a2= -0.00089076;
  b2= 0.0231597;  
  
  G4double ec,ecsq,xnulam,etest(0.),ra(0.),a,w,c,signor(1.),signor2,sig; 
  G4double b,ecut,cut,ecut2,geom,elab;
  
  
  G4double	flow = 1.e-18;
  G4double       spill= 1.e+18; 
  
  
  
  if (ResidualA <= 60.)  signor = 0.92;
  else if (ResidualA < 100.) signor = 0.8 + ResidualA*0.002;
  
  
  ec = 1.44 * theZ * ResidualZ / (1.5*ResidualAthrd+ra);
  ecsq = ec * ec;
  p = p0 + p1/ec + p2/ecsq;
  landa = landa0*ResidualA + landa1;
  a = std::pow(ResidualA,mu1);
  mu = mu0 * a;
  nu = a* (nu0+nu1*ec+nu2*ecsq);
  
  c =std::min(3.15,ec*0.5);
  w = 0.7 * c / 3.15; 
  
  xnulam = nu / landa;
  if (xnulam > spill) xnulam=0.;
  if (xnulam >= flow) etest =std::sqrt(xnulam) + 7.;
  
  a = -2.*p*ec + landa - nu/ecsq;
  b = p*ecsq + mu + 2.*nu/ec;
  ecut = 0.;
  cut = a*a - 4.*p*b;
  if (cut > 0.) ecut = std::sqrt(cut);
  ecut = (ecut-a) / (p+p);
  ecut2 = ecut;
  if (cut < 0.) ecut2 = ecut - 2.;
  elab = K * FragmentA / ResidualA;
  sig = 0.;
  if (elab <= ec) { //start for E<Ec 
    if (elab > ecut2)  sig = (p*elab*elab+a*elab+b) * signor;
    
    signor2 = (ec-elab-c) / w;
    signor2 = 1. + std::exp(signor2);
    sig = sig / signor2;
    // signor2 is empirical global corr'n at low elab for protons in PRECO, not enough for light targets
    //  refinement for proton cross section
    if (ResidualZ<=26)
      sig = sig*std::exp(-(a2*ResidualZ + b2)*(elab-(afit*ResidualZ+bfit)*ec)*(elab-(afit*ResidualZ+bfit)*ec)); 
    
  }              //end for E<Ec
  else{           //start for  E>Ec
    sig = (landa*elab+mu+nu/elab) * signor;
    
    //  refinement for proton cross section
    if ( ResidualZ<=26 && elab <=(afit*ResidualZ+bfit)*ec) 
      sig = sig*std::exp(-(a2*ResidualZ + b2)*(elab-(afit*ResidualZ+bfit)*ec)*(elab-(afit*ResidualZ+bfit)*ec));   
    //
    
    geom = 0.;
    
    if (xnulam < flow || elab < etest) 
      {
        if (sig <0.0) {sig=0.0;}
        return sig;
      }
    geom = std::sqrt(theA*K);
    geom = 1.23*ResidualAthrd + ra + 4.573/geom;
    geom = 31.416 * geom * geom;
    sig = std::max(geom,sig);
    
  }   //end for E>Ec

  return sig;}



//   ************************** end of cross sections ******************************* 



