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
// $Id: G4GammaConversionToMuons.cc 106961 2017-10-31 08:36:29Z gcosmo $
//
//         ------------ G4GammaConversionToMuons physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, April 2002
//
//
// 07-08-02: missprint in OR condition in DoIt : f1<0 || f1>f1_max ..etc ...
// 25-10-04: migrade to new interfaces of ParticleChange (vi)
// ---------------------------------------------------------------------------

#include "G4GammaConversionToMuons.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4EmProcessSubType.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

using namespace std;

static const G4double sqrte=sqrt(exp(1.));
static const G4double PowSat=-0.88;

G4GammaConversionToMuons::G4GammaConversionToMuons(const G4String& processName,
                                                   G4ProcessType type)
  : G4VDiscreteProcess (processName, type),
    Mmuon(G4MuonPlus::MuonPlus()->GetPDGMass()),
    Rc(elm_coupling/Mmuon),
    LowestEnergyLimit (4.*Mmuon), // 4*Mmuon
    HighestEnergyLimit(1e21*eV), // ok to 1e21eV=1e12GeV, then LPM suppression
    CrossSecFactor(1.)
{ 
  SetProcessSubType(fGammaConversionToMuMu);
  MeanFreePath = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4GammaConversionToMuons::~G4GammaConversionToMuons() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4GammaConversionToMuons::IsApplicable(
                                        const G4ParticleDefinition& particle)
{
   return ( &particle == G4Gamma::Gamma() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversionToMuons::BuildPhysicsTable(const G4ParticleDefinition&)
// Build cross section and mean free path tables
{  //here no tables, just calling PrintInfoDefinition
   PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversionToMuons::GetMeanFreePath(const G4Track& aTrack,
                                                   G4double, G4ForceCondition*)

// returns the photon mean free path in GEANT4 internal units
// (MeanFreePath is a private member of the class)

{
   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
   G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
   const G4Material* aMaterial = aTrack.GetMaterial();

   MeanFreePath = (GammaEnergy <= LowestEnergyLimit) 
     ? DBL_MAX : ComputeMeanFreePath(GammaEnergy,aMaterial);

   return MeanFreePath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4GammaConversionToMuons::ComputeMeanFreePath(G4double GammaEnergy,
                                              const G4Material* aMaterial)

// computes and returns the photon mean free path in GEANT4 internal units
{
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

  G4double SIGMA = 0.0;

  for ( size_t i=0 ; i < aMaterial->GetNumberOfElements(); ++i)
  {
    SIGMA += NbOfAtomsPerVolume[i] *
      ComputeCrossSectionPerAtom(GammaEnergy,
                                 (*theElementVector)[i]->GetZasInt());
  }
  return (SIGMA > 0.0) ? 1./SIGMA : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversionToMuons::GetCrossSectionPerAtom(
                                   const G4DynamicParticle* aDynamicGamma,
                                   const G4Element* anElement)

// gives the total cross section per atom in GEANT4 internal units
{
   return ComputeCrossSectionPerAtom(aDynamicGamma->GetKineticEnergy(),
                                     anElement->GetZasInt());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4GammaConversionToMuons::ComputeCrossSectionPerAtom(
                         G4double Egam, G4int Z)
			 
// Calculates the microscopic cross section in GEANT4 internal units.
// Total cross section parametrisation from H.Burkhardt
// It gives a good description at any energy (from 0 to 10**21 eV)
{ 
  if(Egam <= LowestEnergyLimit) return 0.0; // below threshold return 0

  G4double CrossSection = 0.0;
  G4NistManager* nist = G4NistManager::Instance();

  G4double PowThres,Ecor,B,Dn,Zthird,Winfty,WMedAppr,
    Wsatur,sigfac;
  
  if(Z==1) // special case of Hydrogen
    { B=202.4;
      Dn=1.49;
    }
  else
    { B=183.;
      Dn=1.54*nist->GetA27(Z);
    }
  Zthird=1./nist->GetZ13(Z); // Z**(-1/3)
  Winfty=B*Zthird*Mmuon/(Dn*electron_mass_c2);
  WMedAppr=1./(4.*Dn*sqrte*Mmuon);
  Wsatur=Winfty/WMedAppr;
  sigfac=4.*fine_structure_const*Z*Z*Rc*Rc;
  PowThres=1.479+0.00799*Dn;
  Ecor=-18.+4347./(B*Zthird);
  
  G4double CorFuc=1.+.04*G4Log(1.+Ecor/Egam);
  //G4double Eg=pow(1.-4.*Mmuon/Egam,PowThres)*pow( pow(Wsatur,PowSat)+
  //            pow(Egam,PowSat),1./PowSat); // threshold and saturation
  G4double Eg=G4Exp(G4Log(1.-4.*Mmuon/Egam)*PowThres)*
    G4Exp(G4Log( G4Exp(G4Log(Wsatur)*PowSat)+G4Exp(G4Log(Egam)*PowSat))/PowSat);
  CrossSection=7./9.*sigfac*G4Log(1.+WMedAppr*CorFuc*Eg);
  CrossSection*=CrossSecFactor; // increase the CrossSection by  (by default 1)
  return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4GammaConversionToMuons::SetCrossSecFactor(G4double fac)
// Set the factor to artificially increase the cross section
{ 
  CrossSecFactor=fac;
  G4cout << "The cross section for GammaConversionToMuons is artificially "
         << "increased by the CrossSecFactor=" << CrossSecFactor << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VParticleChange* G4GammaConversionToMuons::PostStepDoIt(
                                                        const G4Track& aTrack,
                                                        const G4Step&  aStep)
//
// generation of gamma->mu+mu-
//
{
  aParticleChange.Initialize(aTrack);
  const G4Material* aMaterial = aTrack.GetMaterial();

  // current Gamma energy and direction, return if energy too low
  const G4DynamicParticle *aDynamicGamma = aTrack.GetDynamicParticle();
  G4double Egam = aDynamicGamma->GetKineticEnergy();
  if (Egam <= LowestEnergyLimit) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();

  // select randomly one element constituting the material
  const G4Element* anElement = SelectRandomAtom(aDynamicGamma, aMaterial);
  G4int Z = anElement->GetZasInt();
  G4NistManager* nist = G4NistManager::Instance();

  G4double A027 = nist->GetA27(Z);

  /*
  G4double B,Dn;
  if(Z==1) // special case of Hydrogen
    { B=202.4;
      Dn=1.49;
    }
  else
    { B=183.;
      Dn=1.54*A027;
    }
  G4double Zthird=1./nist->GetZ13(Z); // Z**(-1/3)
  G4double Winfty=B*Zthird*Mmuon/(Dn*electron_mass_c2);

  G4double C1Num=0.138*A027;
  G4double C1Num2=C1Num*C1Num;
  G4double C2Term2=electron_mass_c2/(183.*Zthird*Mmuon);

  G4double GammaMuonInv=Mmuon/Egam;
  G4double sqrtx=sqrt(.25-GammaMuonInv);
  G4double xmax=.5+sqrtx;
  G4double xmin=.5-sqrtx;

  // generate xPlus according to the differential cross section by rejection
  /*  G4double Ds2=(Dn*sqrte-2.);
      G4double sBZ=sqrte*B*Zthird/electron_mass_c2;
      G4double LogWmaxInv=1./G4Log(Winfty*(1.+2.*Ds2*GammaMuonInv)
      /(1.+2.*sBZ*Mmuon*GammaMuonInv));*/
  G4double xPlus,xMinus,xPM;
  G4int nn = 0;
  const G4int nmax = 1000;
  const G4double maxWeight = 1.5; // NOT SURE WHAT THE PARAMETRICS IS!
  // THe studies I have done show that the max weight over 10000 samples can go up to ~1.3 for TeV beam energies; tends to be robustly below 1 for 1-10 GeV beam energies.  In the limit of rho<<u, calculations justify a max weight of 1.  Larger weights can occur for rho~u. 
  // larger weights also seem to occur slightly more frequently with Tsai ME

  G4bool goodEvent = false;
  const G4bool debug = false;

  // NT: This code samples all the parameters: xPlus, t, psi, rho, and then checks a full matrix element
  // (for now, doesn't actually check the ME, just samples and evaluates ME.  I need to pick a reasonable maximum value for the rejection sampling.
  // For now: no form factor at all!



  G4double thetaPlus=0.,thetaMinus=0.,phiHalf=0.; // final angular variables

  // scale of q/mMu at which form factor sets in [PRM(6.115),muonNotes(10)]
  G4double qFormFactor = 1/(0.20*A027); 
  G4double aRangeNorm = log((pow(qFormFactor,4.)+pow(2.*GammaMuonInv,4.))/pow(2.*GammaMuonInv,4.));  // scale for the log(rho) integral.  NOT the best definition (see discussion surrounding (8) and (9) in muonNotes), but within 10-25% for Tungsten at energies of interest to us.  Also worth noting that this normalization factor is only for convenience, so that we can reduce the photon-energy-dependence of the maxweight. 

  do
  { 
    ++nn;

    G4double xPRange=xmax-xmin;

    // Sample xPlus from xmin to xmax
    xPlus=xmin+G4UniformRand()*(xmax-xmin);

    xMinus=1.-xPlus;
    xPM=xPlus*xMinus;
    //    G4double del=Mmuon*Mmuon/(2.*Egam*xPM);
    G4double delOverMmuon=Mmuon/(2.*Egam*xPM);

    /*    W=Winfty*(1.+Ds2*del/Mmuon)/(1.+sBZ*del);
    if(W<=1. || nn > nmax) { break; } // to avoid negative cross section at xmin
    G4double xxp=1.-4./3.*xPM; // the main xPlus dependence
    result=xxp*G4Log(W)*LogWmaxInv;
    if(result>1.) {
      G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
             << " in dSigxPlusGen, result=" << result << " > 1" << G4endl;
    }
    ++nn;
    if(nn >= nmax) { break; }
    */

    // now generate the angular variables via the auxilary variables t,psi,rho

  G4double a3 = (GammaMuonInv/(2.*xPM));
  G4double a33 = a3*a3;
  G4double f1;
  G4double b1  = 1./(4.*C1Num2);
  G4double b3  = b1*b1*b1;
  G4double a21 = a33 + b1;
  
  G4double f1_max=-(1.-xPM)*(2.*b1+(a21+a33)*G4Log(a33/a21))/(2*b3);  

  G4double thetaPlus,thetaMinus,phiHalf; // final angular variables
  nn = 0;
  // t, psi, rho generation start  (while angle < pi)
  do {
    //generate t by the rejection method
    do { 
      ++nn;
      t=G4UniformRand();
      G4double a34=a33/(t*t);
      G4double a22 = a34 + b1;
      if(std::abs(b1)<0.0001*a34) 
	// special case of a34=a22 because of logarithm accuracy
	{
	  f1=(1.-2.*xPM+4.*xPM*t*(1.-t))/(12.*a34*a34*a34*a34);
	}
      else
	{
	  f1=-(1.-2.*xPM+4.*xPM*t*(1.-t))*(2.*b1+(a22+a34)*G4Log(a34/a22))/(2*b3);      
	}
      if(f1<0.0 || f1> f1_max) // should never happend
	{
	  G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
		 << "outside allowed range f1=" << f1 
		 << " is set to zero, a34 = "<< a34 << " a22 = "<<a22<<"."
		 << G4endl;
	  f1 = 0.0;
	}
      if(nn > nmax) { break; }
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko  
    } while ( G4UniformRand()*f1_max > f1);
    // generate psi by the rejection method
    G4double f2_max=1.-2.*xPM*(1.-4.*t*(1.-t));
    // long version
    G4double f2;
    do { 
      ++nn;
      psi=twopi*G4UniformRand();
      f2=1.-2.*xPM+4.*xPM*t*(1.-t)*(1.+cos(2.*psi));
      if(f2<0 || f2> f2_max) // should never happend
	{
	  G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
		 << "outside allowed range f2=" << f2 << " is set to zero"
		 << G4endl;
          f2 = 0.0;
	}
      if(nn >= nmax) { break; }
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while ( G4UniformRand()*f2_max > f2);

    // generate rho by direct transformation
    G4double C2Term1=GammaMuonInv/(2.*xPM*t);
    G4double C22 = C2Term1*C2Term1+C2Term2*C2Term2;
    G4double C2=4.*C22*C22/sqrt(xPM);
    G4double rhomax=(1./t-1.)*1.9/A027;
    G4double beta=G4Log( (C2+rhomax*rhomax*rhomax*rhomax)/C2 );
    rho=G4Exp(G4Log(C2 *( G4Exp(beta*G4UniformRand())-1. ))*0.25);

    G4double rhoMaxForTheta = 2.*u/fabs(cos(psi));
    G4double rhoMaxForPhi = pi*u/fabs(sin(psi));

    G4double rhoMax=fmin(rhoMaxForTheta,rhoMaxForPhi);
    G4double c2 = delOverMmuon/t;

    // Sample a=log(rho^4+c2^4) from rho=0 to rhoMax
    G4double amin=log(pow(c2,4.));
    G4double amax=log(pow(rhoMax,4.)+pow(c2,4.));
    G4double aRange=amax-amin;

    G4double a=amin+G4UniformRand()*aRange;
    G4double rho=pow(G4Exp(a)-pow(c2,4.),0.25);


    //now get the kinematical variables
    G4double xiHalf=0.5*rho*cos(psi);
    G4double phi=rho/u*sin(psi);

    G4double uPlus = u+xiHalf;
    G4double uMinus = u-xiHalf;

    // Only do these calculations if kinematics is sensible...easier to do the u upper limit as a rejection than by changing the rho range
    if(uPlus > xPlus*pi/GammaMuonInv || uMinus > xMinus*pi/GammaMuonInv)
      {
	//	thetaPlus =GammaMuonInv*uPlus/xPlus;
	if(debug) {
	  G4cerr << " [out of range] xP " << xPlus 
		 << " t " << t 
		 << " psi " << psi
		 << " rho " << rho
		 << " uPlus " << uPlus
		 << " uMinus " << uMinus
		 << " thetaPlus " << GammaMuonInv*uPlus/xPlus
		 << " thetaMinus " << GammaMuonInv*uMinus/xMinus << G4endl;
	}
      }
    else
      {
	// All kinematic quantities here are in units of mMuon except where GeV appears in variable names
	G4double gamma0=1/GammaMuonInv;

	G4double gammaPlus = gamma0*xPlus;
	thetaPlus = GammaMuonInv*uPlus/xPlus;
	G4double betaPlus = sqrt(1.-1./(gammaPlus*gammaPlus));
	G4double kdotpPlus = gamma0 * gammaPlus * (1. - betaPlus * cos(thetaPlus));
	G4double pPerpPlus = betaPlus * gammaPlus * sin(thetaPlus);

	G4double gammaMinus = gamma0*xMinus;
	thetaMinus = GammaMuonInv*uMinus/xMinus;
	G4double betaMinus = sqrt(1.-1./(gammaMinus*gammaMinus));
	G4double kdotpMinus = gamma0 * gammaMinus * (1. - betaMinus * cos(thetaMinus));
	G4double pPerpMinus = betaMinus * gammaMinus * sin(thetaMinus);

	// POSITIVE Squared momentum transfer, in units of muon mass and GeV
	//	G4double q2par = pow(delOverMmuon*(1+xMinus*uPlus*uPlus + xPlus*uMinus*uMinus),2.);
	//	G4double q2perp = pow(uPlus - uMinus,2.) + 2 * uPlus * uMinus * (1.-cos(phi));
	G4double qpar = gamma0 - betaMinus * gammaMinus * cos(thetaMinus) 
			     - betaPlus * gammaPlus * cos(thetaPlus);
	G4double q2par = qpar*qpar;
	G4double q2perp = pow(pPerpPlus - pPerpMinus,2.) + 2 * pPerpPlus * pPerpMinus * (1.-cos(phi));

	G4double q2parApprox = pow(delOverMmuon/t,2.);
	G4double q2perpApprox = rho*rho;

	G4double q2 = q2par + q2perp;
	//    G4double q2GeV = q2 * Mmuon * Mmuon; 
    
	// Now evaluate the matrix element, with appropriate measure factors

	G4double integratedMeasure = xPRange * aRange / aRangeNorm;  
	G4double Jacobian_a_Vs_Rho = G4Exp(a) / (rho*rho*rho);
	G4double Jacobian_tRhoPsi_VS_uPuMPhi = rho/(t*(1.-t));
	G4double fullJacobian = Jacobian_a_Vs_Rho * Jacobian_tRhoPsi_VS_uPuMPhi;

	
	//  OLD: Matrix Element from G4 physics reference manual

	G4double dPlus = 1+uPlus*uPlus;
	G4double dMinus = 1+uMinus*uMinus;

	G4double curlyBrace5_50_A = (uPlus*uPlus + uMinus*uMinus)/(dPlus * dMinus);
	G4double curlyBrace5_50_B = -2.*xPlus*xMinus * 
	  (uPlus*uPlus/(dPlus*dPlus) + uMinus*uMinus/(dMinus*dMinus));
	G4double curlyBrace5_50_C = -2.*uPlus * uMinus* (1-2.*xPlus*xMinus)*cos(phi)/(dPlus*dMinus);
	G4double curlyBrace5_50 = curlyBrace5_50_A + curlyBrace5_50_B + curlyBrace5_50_C;
	G4double matrixEl5_50 = 1./(q2*q2) * uPlus * uMinus * curlyBrace5_50; 
	
	// This matrix element comes from Tsai 1974 paper and should be exact up to 1/m_{Nucleus} corrections (which aren't hard to add) and inelastic terms (also possible to add if we know the appropriate weight and how to sample it...)
	G4double term1Plus  = -(-q2+kdotpMinus)/(kdotpPlus);
	G4double term1Minus = -(-q2+kdotpPlus)/(kdotpMinus);
	G4double term2Plus  = (-q2/2. +2.* gammaMinus * gammaMinus)/(kdotpPlus * kdotpPlus);
	G4double term2Minus = (-q2/2. +2.* gammaPlus * gammaPlus)/(kdotpMinus * kdotpMinus);
	G4double mixedTerm = 2.*((1+q2/2.)*(-2.*gammaMinus * gammaPlus - q2/2.)+q2*gamma0*gamma0/2.)/(kdotpMinus * kdotpPlus);

	G4double sumTerms = term1Plus + term1Minus + term2Plus + term2Minus + mixedTerm;
	G4double matrixElTsai = 1./(4.*q2*q2) * sin(thetaPlus) * sin(thetaMinus) * sumTerms;

	
	G4double formFactorSqrt = 1./(1.+q2/(qFormFactor*qFormFactor));
	G4double formFactor = formFactorSqrt*formFactorSqrt;
	
	G4double theWeight = integratedMeasure * fullJacobian * matrixElTsai * formFactor;
	
	if(theWeight < 0 || theWeight > maxWeight)
	  {
	    G4cerr << " **** WARNING in G4GammaConversionToMuons.cc: weight " <<theWeight 
		   << " out of bounds 0 to " << maxWeight << G4endl;
	  }
	else if(theWeight > 1.0)
	  {
	    G4cerr << " **** WARNING in G4GammaConversionToMuons.cc: weight " <<theWeight 
		   << " is > 1! (This is just informational -- weights up to " << maxWeight << " are handled correctly)"<< G4endl;
	  }
	

	if(debug || theWeight < 0 || theWeight > maxWeight || theWeight > 1.0 ) {
	  G4cerr << " xP " << xPlus 
		 << " t " << t 
		 << " psi " << psi
		 << " rho " << rho
		 << " q2 " << q2
		 << " 1+up2 " << 1+ uPlus*uPlus
		 << " 2k.pxp " << 2.*kdotpPlus*xPlus
		 << " 1+um2 " << 1+ uMinus*uMinus
		 << " 2k.pxm " << 2.*kdotpMinus*xMinus
		


		 << "thetaP " << thetaPlus
		 << " thetaM " << thetaMinus
		 << " q2Approx " << q2perpApprox + q2parApprox
	    //		 << " curly " << curlyBrace5_50
		 << " ME5.50 " << matrixEl5_50
		 << " METsai " << matrixElTsai
		 << " terms " << term1Plus << " " << term1Minus << " " << term2Plus << " " << term2Minus << " " << mixedTerm 
		 << " measurr " << integratedMeasure 
		 << " ff " << formFactor
		 << " weight " << theWeight
		 << G4endl;
	}
	
	if(G4UniformRand() * maxWeight < theWeight) 
	  goodEvent = true;
	
	phiHalf=0.5*phi;
	//	thetaPlus =GammaMuonInv*(u+xiHalf)/xPlus;
	//	thetaMinus=GammaMuonInv*(u-xiHalf)/xMinus;
      }
    

  } while (nn<nmax && !goodEvent) ;

  if(!goodEvent) {
    G4cerr << "**** WARNING in G4GammaConversionToMuons.cc: failed to select a ME-weighted event after " << nmax << " tries.  Using unphysical event weight." << G4endl;
  }

  // now construct the vectors
  // azimuthal symmetry, take phi0 at random between 0 and 2 pi
  G4double phi0=twopi*G4UniformRand(); 
  G4double EPlus=xPlus*Egam;
  G4double EMinus=xMinus*Egam;

  // mu+ mu- directions for gamma in z-direction
  G4ThreeVector MuPlusDirection  ( sin(thetaPlus) *cos(phi0+phiHalf),
                   sin(thetaPlus)  *sin(phi0+phiHalf), cos(thetaPlus) );
  G4ThreeVector MuMinusDirection (-sin(thetaMinus)*cos(phi0-phiHalf),
                  -sin(thetaMinus) *sin(phi0-phiHalf), cos(thetaMinus) );
  // rotate to actual gamma direction
  MuPlusDirection.rotateUz(GammaDirection);
  MuMinusDirection.rotateUz(GammaDirection);
  aParticleChange.SetNumberOfSecondaries(2);
  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1= new G4DynamicParticle(
                           G4MuonPlus::MuonPlus(),MuPlusDirection,EPlus-Mmuon);
  aParticleChange.AddSecondary(aParticle1);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2= new G4DynamicParticle(
                       G4MuonMinus::MuonMinus(),MuMinusDirection,EMinus-Mmuon);
  aParticleChange.AddSecondary(aParticle2);
  //
  // Kill the incident photon
  //
  aParticleChange.ProposeMomentumDirection( 0., 0., 0. ) ;
  aParticleChange.ProposeEnergy( 0. ) ;
  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
  //  Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element* G4GammaConversionToMuons::SelectRandomAtom(
		  const G4DynamicParticle* aDynamicGamma,
		  const G4Material* aMaterial)
{
  // select randomly 1 element within the material, invoked by PostStepDoIt

  const G4int NumberOfElements            = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4Element* elm = (*theElementVector)[0];

  if (NumberOfElements > 1) { 
    const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

    G4double PartialSumSigma = 0.;
    G4double rval = G4UniformRand()/MeanFreePath;

    for (G4int i=0; i<NumberOfElements; ++i)
    { 
      elm = (*theElementVector)[i];
      PartialSumSigma += NbOfAtomsPerVolume[i]
	*GetCrossSectionPerAtom(aDynamicGamma, elm);
      if (rval <= PartialSumSigma) { break; }
    }
  }
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4GammaConversionToMuons::PrintInfoDefinition()
{
  G4String comments ="gamma->mu+mu- Bethe Heitler process, SubType= ";
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << GetProcessSubType() << G4endl;
  G4cout << "        good cross section parametrization from "
         << G4BestUnit(LowestEnergyLimit,"Energy")
         << " to " << HighestEnergyLimit/GeV << " GeV for all Z." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
