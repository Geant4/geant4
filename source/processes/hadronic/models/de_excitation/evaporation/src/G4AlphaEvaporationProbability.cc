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
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option

#include "G4AlphaEvaporationProbability.hh"

G4AlphaEvaporationProbability::G4AlphaEvaporationProbability() :
    G4EvaporationProbability(4,2,1,&theCoulombBarrier) // A,Z,Gamma,&theCoumlombBarrier
{
	
}


G4AlphaEvaporationProbability::G4AlphaEvaporationProbability(const G4AlphaEvaporationProbability &): G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4AlphaEvaporationProbability & G4AlphaEvaporationProbability::
operator=(const G4AlphaEvaporationProbability &) 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4AlphaEvaporationProbability::operator==(const G4AlphaEvaporationProbability &) const
{
    return false;
}

G4bool G4AlphaEvaporationProbability::operator!=(const G4AlphaEvaporationProbability &) const
{
    return true;
}


 G4double G4AlphaEvaporationProbability::CalcAlphaParam(const G4Fragment & fragment)  
  { return 1.0 + CCoeficient(static_cast<G4double>(fragment.GetZ()-GetZ()));}
	
 G4double G4AlphaEvaporationProbability::CalcBetaParam(const G4Fragment &) 
  { return 0.0; }

G4double G4AlphaEvaporationProbability::CCoeficient(const G4double aZ) 
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
    G4double C = 0.0;
	
	
    if (aZ <= 30) {
	C = 0.10;
    } else if (aZ <= 50) {
	C = 0.1 + -((aZ-50.)/20.)*0.02;
    } else if (aZ < 70) {
	C = 0.08 + -((aZ-70.)/20.)*0.02;
    } else {
	C = 0.06;
    }
    return C;
}

///////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
 G4double G4AlphaEvaporationProbability::CrossSection(const  G4Fragment & fragment, const  G4double K)
{
  theA=GetA();
  theZ=GetZ();
  ResidualA=fragment.GetA()-theA;
  ResidualZ=fragment.GetZ()-theZ; 
  
  ResidualAthrd=std::pow(ResidualA,0.33333);
  FragmentA=fragment.GetA();
  FragmentAthrd=std::pow(FragmentA,0.33333);
  
  
  if (OPTxs==0) {std::ostringstream errOs;
    errOs << "We should'n be here (OPT =0) at evaporation cross section calculation (Alpha's)!!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;}
  
  if( OPTxs==1 || OPTxs==2) return G4AlphaEvaporationProbability::GetOpt12( K);
  else if (OPTxs==3 || OPTxs==4)  return G4AlphaEvaporationProbability::GetOpt34( K);
  else{
    std::ostringstream errOs;
    errOs << "BAD Alpha CROSS SECTION OPTION AT EVAPORATION!!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
}


//
//********************* OPT=1,2 : Chatterjee's cross section ************************ 
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4AlphaEvaporationProbability::GetOpt12(const  G4double K)
{
// c     ** alpha from huizenga and igo
  G4double Kc=K;
  
  // JMQ xsec is set constat above limit of validity
  if (K>50) Kc=50;
  
  G4double landa ,mu ,nu ,p , Ec,q,r,ji,xs;
  
  G4double     p0 = 10.95;
  G4double     p1 = -85.2;
  G4double     p2 = 1146.;
  G4double     landa0 = 0.0643;
  G4double     landa1 = -13.96;
  G4double     mu0 = 781.2;
  G4double     mu1 = 0.29;
  G4double     nu0 = -304.7;
  G4double     nu1 = -470.0;
  G4double     nu2 = -8.580;   
  G4double     delta=1.2;          

 
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
G4double G4AlphaEvaporationProbability::GetOpt34(const  G4double K)
// c     ** alpha from huizenga and igo
{
  
  G4double landa, mu, nu, p , signor(1.),sig;
  G4double ec,ecsq,xnulam,etest(0.),a; 
  G4double b,ecut,cut,ecut2,geom,elab;


  G4double     flow = 1.e-18;
  G4double     spill= 1.e+18;

  G4double       p0 = 10.95;
  G4double     p1 = -85.2;
  G4double     p2 = 1146.;
  G4double     landa0 = 0.0643;
  G4double     landa1 = -13.96;
  G4double     mu0 = 781.2;
  G4double     mu1 = 0.29;
  G4double     nu0 = -304.7;
  G4double     nu1 = -470.0;
  G4double     nu2 = -8.580;        
  
  G4double      ra=1.20;
  
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
