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
//J.M.Quesada (August 08). New source file
//
// Modif (21 August 2008) by J. M. Quesada for external choice of inverse 
// cross section option
 
#include "G4PreCompoundTriton.hh"


  G4ReactionProduct * G4PreCompoundTriton::GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct =
      new G4ReactionProduct(G4Triton::TritonDefinition());
    theReactionProduct->SetMomentum(GetMomentum().vect());
    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef PRECOMPOUND_TEST
    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
    return theReactionProduct;
  }   

   G4double G4PreCompoundTriton::FactorialFactor(const G4double N, const G4double P)
  {
      return 
      (N-3.0)*(P-2.0)*(
		       (((N-2.0)*(P-1.0))/2.0) *(
						 (((N-1.0)*P)/3.0) 
						 )
		       );
  }
  
   G4double G4PreCompoundTriton::CoalescenceFactor(const G4double A)
  {
     return 243.0/(A*A);
  }    


   G4double G4PreCompoundTriton::GetRj(const G4int NumberParticles, const G4int NumberCharged)
  {
    G4double rj = 0.0;
    G4double denominator = NumberParticles*(NumberParticles-1)*(NumberParticles-2);
    if(NumberCharged >= 1 && (NumberParticles-NumberCharged) >= 2) rj = 3.0*static_cast<G4double>(NumberCharged*(NumberParticles-NumberCharged)*(NumberParticles-NumberCharged-1))/static_cast<G4double>(denominator); 
    return rj;
  }




////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
 G4double G4PreCompoundTriton::CrossSection(const  G4double K)
{
  ResidualA=GetRestA();
  ResidualZ=GetRestZ(); 
  theA=GetA();
  theZ=GetZ();
  ResidualAthrd=std::pow(ResidualA,0.33333);
  FragmentA=GetA()+GetRestA();
  FragmentAthrd=std::pow(FragmentA,0.33333);

  if (OPTxs==0) return GetOpt0( K);
  else if( OPTxs==1 || OPTxs==2) return GetOpt12( K);
  else if (OPTxs==3 || OPTxs==4)  return GetOpt34( K);
  else{
    std::ostringstream errOs;
    errOs << "BAD TRITON CROSS SECTION OPTION !!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
}

// *********************** OPT=0 : Dostrovski's cross section  *****************************

G4double G4PreCompoundTriton::GetOpt0(const  G4double K)
{
  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  //with the rest of the options
  return 1.e+25*pi*(r0*ResidualAthrd)*(r0*ResidualAthrd)*GetAlpha()*(1.+GetBeta()/K);
}
//
//---------
//
  G4double G4PreCompoundTriton::GetAlpha()
  {
    G4double C = 0.0;
    G4double aZ = GetZ() + GetRestZ();
    if (aZ >= 70) 
      {
	C = 0.10;
      } 
    else 
      {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375; 
      }
 
    return 1.0 + C/3.0;
  }
//
//-------------
//
   G4double G4PreCompoundTriton::GetBeta() 
  {
      return -GetCoulombBarrier();
  }
//
//********************* OPT=1,2 : Chatterjee's cross section ************************ 
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundTriton::GetOpt12(const  G4double K)
{

  G4double Kc=K;

  // JMQ xsec is set constat above limit of validity
  if (K>50) Kc=50;

  G4double landa ,mu ,nu ,p , Ec,q,r,ji,xs;
 
  G4double    p0 = -11.04;
  G4double    p1 = 619.1;
  G4double    p2 = -2147.;
  G4double    landa0 = -0.0426;
  G4double    landa1 = -10.33;
  G4double    mu0 = 601.9;
  G4double    mu1 = 0.37;
  G4double    nu0 = 583.0;
  G4double    nu1 = -546.2;
  G4double    nu2 = 1.718;  
  G4double    delta=1.2;            

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

// *********** OPT=3,4 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundTriton::GetOpt34(const  G4double K)
//     ** t from o.m. of hafele, flynn et al
{

  G4double landa, mu, nu, p , signor(1.),sig;
G4double ec,ecsq,xnulam,etest(0.),a; 
G4double b,ecut,cut,ecut2,geom,elab;


 G4double     flow = 1.e-18;
 G4double     spill= 1.e+18;


 G4double     p0 = -21.45;
 G4double     p1 = 484.7;
 G4double     p2 = -1608.;
 G4double     landa0 = 0.0186;
 G4double     landa1 = -8.90;
 G4double     mu0 = 686.3;
 G4double     mu1 = 0.325;
 G4double     nu0 = 368.9;
 G4double     nu1 = -522.2;
 G4double     nu2 = -4.998;  
  
 G4double      ra=0.80;
        
 ec = 1.44 * theZ * ResidualZ / (1.5*ResidualAthrd+ra);
 ecsq = ec * ec;
 p = p0 + p1/ec + p2/ecsq;
 landa = landa0*ResidualA + landa1;
 a = std::pow(ResidualA,mu1);
 mu = mu0 * a;
 nu = a* (nu0+nu1*ec+nu2*ecsq);  
 xnulam = nu / landa;
 if (xnulam > spill) xnulam=0.;
 if (xnulam >= flow) etest = 1.2 *std::sqrt(xnulam);
 
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
 }           //end for E<Ec
 else {           //start for E>Ec
   sig = (landa*elab+mu+nu/elab) * signor;
   geom = 0.;
   if (xnulam < flow || elab < etest) return sig;
   geom = std::sqrt(theA*K);
   geom = 1.23*ResidualAthrd + ra + 4.573/geom;
   geom = 31.416 * geom * geom;
   sig = std::max(geom,sig);
 }           //end for E>Ec
 return sig;

}

//   ************************** end of cross sections ******************************* 




