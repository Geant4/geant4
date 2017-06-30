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
// $Id: G4GEMChannelVI.cc 98577 2016-07-25 13:05:12Z vnivanch $
//
// GEM de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4GEMChannelVI.hh"
#include "G4VCoulombBarrier.hh"
#include "G4CoulombBarrier.hh"
#include "G4PairingCorrection.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

// 10-Points Gauss-Legendre abcisas and weights
const G4double G4GEMChannelVI::ws[] = {
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
const G4double G4GEMChannelVI::xs[] = {
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

G4GEMChannelVI::G4GEMChannelVI(G4int theA, G4int theZ)
  : A(theA), Z(theZ), levelDensity(0.1)
{ 
  fG4pow = G4Pow::GetInstance(); 
  Z13 = fG4pow->Z13(Z);
  A13 = fG4pow->Z13(A);

  cBarrier = new G4CoulombBarrier(A, Z);
  pairingCorrection = G4PairingCorrection::GetInstance();

  nData = G4NuclearLevelData::GetInstance();
  levelManager = nData->GetLevelManager(Z, A);
  maxLevelE = levelManager->MaxLevelEnergy();

  massGround = G4NucleiProperties::GetNuclearMass(A, Z);

  resA = resZ = fragZ = fragA = nWarn = 0;
  massGround = maxLevelE = Z13 = A13 = massFrag = eCBarrier 
    = resMassGround = resZ13 = resA13 = delta0
    = delta1 = maxExc = maxProb = alphaP = betaP = maxKinEnergy = 0.0;
  coeff = CLHEP::fermi*CLHEP::fermi/(CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);
}

G4GEMChannelVI::~G4GEMChannelVI()
{
  delete cBarrier;
}

void G4GEMChannelVI::Initialise()
{
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  levelDensity = param->GetLevelDensity();
}

G4double G4GEMChannelVI::GetEmissionProbability(G4Fragment* fragment)
{
  fragZ = fragment->GetZ_asInt();
  fragA = fragment->GetA_asInt();
  resZ = fragZ - Z;
  resA = fragA - A;
  G4double prob = 0.0;
  if(resA < A || resA < resZ || resZ < 0 || (resA == A && resZ < Z)) { 
    return prob; 
  }
 
  resMassGround = G4NucleiProperties::GetNuclearMass(resA, resZ);
  G4double exc = fragment->GetExcitationEnergy();
  massFrag = fragment->GetGroundStateMass() + exc;
  delta0 = pairingCorrection->GetPairingCorrection(fragA, fragZ);
  eCBarrier = cBarrier->GetCoulombBarrier(resA, resZ, exc);
  
  maxExc = massFrag - massGround - resMassGround - eCBarrier - delta0;
  if(maxExc < 0.0) { return prob; }

  resZ13 = fG4pow->Z13(resZ);
  resA13 = fG4pow->Z13(resA);
  delta1 = pairingCorrection->GetPairingCorrection(resA, resZ);

  G4double C = 0.0;
  if(resA >= 50) {
    C = -0.10/G4double(A);
  } else if(resZ > 20) {
    C = (0.123482-0.00534691*Z-0.0000610624*(Z*Z)+5.93719*1e-7*(Z*Z*Z)+
	 1.95687*1e-8*(Z*Z*Z*Z))/G4double(A);
  }
  if(0 == Z) {
    alphaP = 0.76+1.93/resA13;
    betaP = (1.66/(resA13*resA13)-0.05)*CLHEP::MeV/alphaP;
  } else {
    alphaP = 1.0 + C;
    betaP = - eCBarrier;
  }

  maxProb = 0.0;
  G4double e0 = maxExc*0.5;

  // e is an excitation of emitted fragment 
  for (G4int i=0; i<NPOINTSGEM; ++i) {
    prob += ws[i]*IntegratedProbability(e0*(xs[i] + 1.0));
  }
  prob *= coeff*e0*e0;

  /*
  G4cout << "G4GEMChannelVI: Z= " << Z << "  A= " << A 
	 << " FragmentZ= " << aZ << " FragmentA= " << anA
	 << " Zres= " << ResidualZ << " Ares= " << ResidualA 
	 << G4endl; 
  */
  
  //G4cout << "Prob= " << prob << G4endl;
  return prob;
}

G4double G4GEMChannelVI::IntegratedProbability(G4double exc)
{
  G4double e0 = (maxExc - exc)*0.5;

  G4double y;
  G4double sum = 0.0;

  for (G4int i=0; i<NPOINTSGEM; ++i) {
    y = ProbabilityDistributionFunction(exc, e0*(xs[i] + 1.0));
    maxProb = std::max(maxProb, y);
    sum += ws[i]*y;
  }
  return sum;
}

G4Fragment* G4GEMChannelVI::EmittedFragment(G4Fragment* theNucleus)
{
  G4double exc;
  G4double resExc;

  static G4double factor = 1.2;
  maxProb *= factor;
 
  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  for(G4int i=0; i<100; ++i) {
    do {
      exc = maxExc*rndm->flat();
      resExc = maxExc*rndm->flat();
    } while (exc + resExc > maxExc);

    G4double prob = ProbabilityDistributionFunction(exc, resExc);
    if(prob > maxProb && nWarn < 10) {
      ++nWarn;
      G4cout << "### G4GEMChannelVI::EmittedFragment WARNING: majoranta "
	     << maxProb << " is exceeded " << prob << "\n"
	     << " fragZ= " << fragZ << " fragA= " << fragA
	     << " Z= " << Z << " A= " << A 
	     << " resZ= " << resZ << " resA= " << resA << "\n"
             << " exc(MeV)= " << exc << " resExc(MeV)= " << resExc 
	     << " maxExc(MeV)= " << maxExc << G4endl;
    }
    if(maxProb*rndm->flat() <= prob) { break; }
  }
  if(exc <= maxLevelE) {
    exc = FindLevel(levelManager, exc, maxExc - resExc);
  }
  if(resA >= nData->GetMinA(resZ) && resA <= nData->GetMaxA(resZ) 
     && resExc < nData->GetMaxLevelEnergy(Z, A)) {
    const G4LevelManager* lman = nData->GetLevelManager(Z, A);
    if(lman) { resExc = FindLevel(lman, resExc, maxExc - exc); }
  }

  G4LorentzVector lv0 = theNucleus->GetMomentum();
  G4double mass1 = massGround + exc;
  G4double mass2 = resMassGround + resExc;

  G4double e1 = 0.5*((massFrag - mass2)*(massFrag + mass2) 
		     + mass1*mass1)/massFrag;

  G4double p1(0.0);
  if(e1 > mass1) {
    p1 = std::sqrt((e1 - mass1)*(e1 + mass1));
  } else {
    e1 = mass1;
  }
  G4ThreeVector v = G4RandomDirection();
  G4LorentzVector lv1(p1*v, e1);

  G4ThreeVector boostVector = lv0.boostVector();
  lv1.boost(boostVector);
  
  G4Fragment* frag = new G4Fragment(A, Z, lv1);

  G4double e2 = massFrag - e1;
  if(e2 < mass2) {
    e2 = mass2;
    p1 = 0.0;
  }
  lv0.set(-v*p1, e2);
  lv0.boost(boostVector);

  theNucleus->SetZandA_asInt(resZ, resA);
  theNucleus->SetMomentum(lv0);

  return frag; 
} 

G4double 
G4GEMChannelVI::ProbabilityDistributionFunction(G4double exc, G4double resExc)
{
  G4double Ux = (2.5 + 150.0/G4double(resA))*CLHEP::MeV;
  G4double Ex = Ux + delta1;
  G4double T  = 1.0/(std::sqrt(levelDensity/Ux) - 1.5/Ux);
  G4double E0 = Ex - T*(G4Log(T) - G4Log(levelDensity)*0.25 
	- 1.25*G4Log(Ux) + 2.0*std::sqrt(levelDensity*Ux));

  G4double UxCN = (2.5 + 150.0/G4double(A))*CLHEP::MeV;
  G4double ExCN = UxCN + delta0;
  G4double TCN  = 1.0/(std::sqrt(levelDensity/UxCN) - 1.5/UxCN);

  G4double mass1 = massGround + exc;
  G4double mass2 = resMassGround + resExc;

  maxKinEnergy = 0.5*((massFrag - mass2)*(massFrag + mass2) 
		      + mass1*mass1)/massFrag - mass1;
  maxKinEnergy = std::max(maxKinEnergy, 0.0);

  G4double Width = 0.0;
  G4double t = maxKinEnergy/T;
  if ( maxKinEnergy < Ex ) {
    Width = (I1(t,t)*T + (betaP+eCBarrier)*I0(t))/G4Exp(E0/T);

  } else {

    G4double tx = Ex/T;
    G4double s0 = 2.0*std::sqrt(levelDensity*(maxKinEnergy-delta0));
    G4double sx = 2.0*std::sqrt(levelDensity*(Ex-delta0));

    // VI: protection against FPE exception
    if(s0 > 350.) { s0 = 350.; }

    G4double expE0T = G4Exp(E0/T);
    G4double exps0  = G4Exp(s0);
    static const G4double sqrt2 = std::sqrt(2.0);

    Width = I1(t,tx)*T/expE0T + I3(s0,sx)*exps0/(sqrt2*levelDensity);

    if (0 == Z) {
      Width += (betaP+eCBarrier)*(I0(tx)/expE0T + 2.0*sqrt2*I2(s0,sx)*exps0);
    }
  }
  Width *= alphaP*massFrag;

  
  //JMQ 190709 fix on Rb and  geometrical cross sections according to 
  //           Furihata's paper (JAERI-Data/Code 2001-105, p6)
  G4double Rb = 0.0;
  if (A > 4) {
    Rb = 1.12*(resA13 + A13) - 0.86*((resA13 + A13)/(resA13*A13))+2.85;
  } else if (A > 1) {
    Rb=1.5*(resA13 + A13);
  } else {
    Rb = 1.5*resA13;
  }

  G4double ild;
  if (exc < ExCN ) {
    G4double E0CN = ExCN - TCN*(G4Log(TCN) - 0.25*G4Log(levelDensity) 
				- 1.25*G4Log(UxCN) 
				+ 2.0*std::sqrt(levelDensity*UxCN));
    ild = G4Exp((exc-E0CN)/TCN)/TCN;
  } else {
    G4double x  = exc - delta0;
    G4double x1 = std::sqrt(levelDensity*x);
    ild = G4Exp(2*x1)/(x*std::sqrt(x1));
  }

  Width *= (Rb*Rb/ild); 
  return Width;
}

G4double G4GEMChannelVI::FindLevel(const G4LevelManager* man,
				   G4double exc, G4double exclim)
{
  size_t idx = man->NearestLowEdgeLevelIndex(exc);
  size_t idxm = man->NumberOfTransitions();
  G4double e1 = man->LevelEnergy(idx);
  if(idx + 1 < idxm) {
    G4double e2 = man->LevelEnergy(idx+1);
    if(e2 <= exclim) {
      G4int s1 = std::abs(man->SpinParity(idx))+1;
      G4int s2 = std::abs(man->SpinParity(idx+1))+1;
      G4double pr = (G4double)s1/(G4double)(s1 + s2);
      pr = (exc - e1 <= e2 - exc) ? 1.0 - (1.0 - pr)*2*(exc - e1)/(e2 - e1) :
	2*pr*(e2 - exc)/(e2 - e1);
      exc = (G4UniformRand() < pr) ? e1 : e2; 
    } else {
      exc = e1;
    }
  } else { exc = e1; }
  return exc;
}

G4double G4GEMChannelVI::I3(G4double s0, G4double sx)
{
  G4double s2 = s0*s0;
  G4double sx2 = sx*sx;
  G4double S = 1.0/std::sqrt(s0);
  G4double S2 = S*S;
  G4double Sx = 1.0/std::sqrt(sx);
  G4double Sx2 = Sx*Sx;
  
  G4double p1 = S *(2.0 + S2 *( 4.0 + S2 *( 13.5 + S2 *( 60.0 + S2 * 325.125 ))));
  G4double p2 = Sx*Sx2 *((s2-sx2) 
			 + Sx2 *((1.5*s2+0.5*sx2) 
				 + Sx2 *((3.75*s2+0.25*sx2) 
					 + Sx2 *((12.875*s2+0.625*sx2) 
						 + Sx2 *((59.0625*s2+0.9375*sx2) 
							 + Sx2 *(324.8*s2+3.28*sx2))))));
  p2 *= G4Exp(sx-s0);
  return p1-p2; 
}

void G4GEMChannelVI::Dump() const
{
}



