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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//
#include "G4GEMProbabilityVI.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4PairingCorrection.hh"
#include "G4NucleiProperties.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"

// 10-Points Gauss-Legendre abcisas and weights
/*
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
*/

G4GEMProbabilityVI::G4GEMProbabilityVI(G4int anA, G4int aZ, const G4LevelManager* p) 
  : G4VEmissionProbability(aZ, anA), lManager(p)
{
  fragA  = fragZ = 0;
  resA13 = U = delta0 = delta1 = a0 = a1 = probmax = alphaP = betaP = 0.0;
  Umax = bCoulomb = 0.0;
  Gamma = 1.0;
  pcoeff = Gamma*pEvapMass*CLHEP::millibarn
    /((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 
  coeff = CLHEP::fermi*CLHEP::fermi/(CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  isExcited = (!lManager || 0.0 == lManager->MaxLevelEnergy()) ? false : true;
  A13 = pG4pow->Z13(theA);

  if(0 == aZ) {
    ResetIntegrator(30, 0.25*CLHEP::MeV, 0.02);
  } else {
    ResetIntegrator(30, 0.5*CLHEP::MeV, 0.03);
  }
}

G4GEMProbabilityVI::~G4GEMProbabilityVI() 
{}

G4double G4GEMProbabilityVI::ComputeTotalProbability(
         const G4Fragment& fragment, G4double CB)
{
  fragA = fragment.GetA_asInt();
  fragZ = fragment.GetZ_asInt();

  bCoulomb = CB;
  U = fragment.GetExcitationEnergy();
  delta0 = pNuclearLevelData->GetPairingCorrection(fragZ,fragA);
  delta1 = pNuclearLevelData->GetPairingCorrection(resZ,resA);
  Umax = pMass - pEvapMass - pResMass - CB;
  if(0.0 >= Umax) { return 0.0; }

  resA13 = pG4pow->Z13(resA);
  a0 = pNuclearLevelData->GetLevelDensity(fragZ,fragA,U);

  G4double C = 0.0;
  G4int Z2 = theZ*theZ;
  G4int Z3 = Z2*theZ;
  G4int Z4 = Z2*Z2;

  if(resA >= 50) {
    C = -0.10/(G4double)theA;
  } else if(resZ > 20) {
    C = (0.123482-0.00534691*theZ-0.0000610624*Z2+5.93719*1e-7*Z3+
	 1.95687*1e-8*Z4)/(G4double)theA;
  }
  if(0 == theZ) {
    alphaP = 0.76+1.93/resA13;
    betaP = (1.66/(resA13*resA13)-0.05)*CLHEP::MeV/alphaP;
  } else {
    alphaP = 1.0 + C;
    betaP = - bCoulomb;
  }
  if(isExcited) {
    pProbability = Integrated2DProbability();

  } else {
    const G4double twoMass = pMass + pMass;
    const G4double evapMass2 = pEvapMass*pEvapMass;
    G4double ekinmax = 
     ((pMass-pResMass)*(pMass+pResMass) + evapMass2)/twoMass - pEvapMass;
    G4double ekinmin = 
      std::max((CB*(twoMass - CB) + evapMass2)/twoMass - pEvapMass,0.0);
    if(ekinmax <= ekinmin) { return 0.0; }
    pProbability = IntegrateProbability(ekinmin, ekinmax, CB);
  }
  /*  
  G4cout << "G4GEMProbabilityVI: Z= " << theZ << " A= " << theA 
	 << " resZ= " << resZ << " resA= " << resA 
	 << " fragZ= " << fragZ << " fragA= " << fragA 
         << " prob= " << pProbability 
	 << "\n   U= " << U << " Umax= " << Umax << " d0= " << delta0 
         << " a0= " << a0 << G4endl;
  */
  return pProbability;
}

G4double G4GEMProbabilityVI::ComputeProbability(G4double ekin, G4double)
{ 
  // abnormal case - should never happens
  if(pMass < pEvapMass + pResMass) { return 0.0; }
    
  const G4double m02   = pMass*pMass;
  const G4double m12   = pEvapMass*pEvapMass;
  const G4double mres  = std::sqrt(m02 + m12 - 2.*pMass*(pEvapMass + ekin));

  G4double excRes = std::max(mres - pResMass, 0.0);
  a1 = pNuclearLevelData->GetLevelDensity(resZ,resA,excRes);
  G4double prob = ProbabilityDistributionFunction(0.0, excRes);

  //G4cout<<"### G4GEMProbabilityVI::ComputeProbability: Ekin(MeV)= "<<ekin 
  //<< " excRes(MeV)= " << excRes << " prob= " << prob << << G4endl;
  return prob;
}

G4Fragment* G4GEMProbabilityVI::SampleEvaporationFragment()
{
  if(isExcited) { return Sample2DDistribution(); } 
  G4double ekin = SampleEnergy();
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*pEvapMass))
                     *G4RandomDirection(), ekin + pEvapMass);
  G4Fragment* evFragment = new G4Fragment(theA, theZ, lv);
  return evFragment;
}

G4double G4GEMProbabilityVI::Integrated2DProbability()
{
  return 0.0;
}

G4double G4GEMProbabilityVI::ProbabilityDistributionFunction(
         G4double exc, G4double resExc)
{
  G4double Ux = (2.5 + 150.0/G4double(resA))*CLHEP::MeV;
  G4double Ex = Ux + delta1;
  G4double T  = 1.0/(std::sqrt(a0/Ux) - 1.5/Ux);
  G4double E0 = Ex - T*(G4Log(T) - G4Log(a0)*0.25 
	- 1.25*G4Log(Ux) + 2.0*std::sqrt(a0*Ux));

  G4double UxCN = (2.5 + 150.0/(G4double)theA)*CLHEP::MeV;
  G4double ExCN = UxCN + delta0;
  G4double TCN  = 1.0/(std::sqrt(a0/UxCN) - 1.5/UxCN);

  G4double mass1 = pEvapMass + exc;
  G4double mass2 = pResMass + resExc;

  G4double maxKinEnergy = std::max(0.5*((pMass - mass2)*(pMass + mass2) 
					+ mass1*mass1)/pMass - mass1, 0.0);

  G4double Width = 0.0;
  G4double t = maxKinEnergy/T;
  if ( maxKinEnergy < Ex ) {
    Width = (I1(t,t)*T + (betaP+bCoulomb)*I0(t))/G4Exp(E0/T);

  } else {

    G4double tx = Ex/T;
    G4double s0 = 2.0*std::sqrt(a0*(maxKinEnergy-delta0));
    G4double sx = 2.0*std::sqrt(a0*(Ex-delta0));

    // VI: protection against FPE exception
    s0 = std::min(s0, 350.);

    G4double expE0T = G4Exp(E0/T);
    G4double exps0  = G4Exp(s0);
    const G4double sqrt2 = std::sqrt(2.0);

    Width = I1(t,tx)*T/expE0T + I3(s0,sx)*exps0/(sqrt2*a0);

    if (0 == theZ) {
      Width += (betaP+bCoulomb)*(I0(tx)/expE0T + 2.0*sqrt2*I2(s0,sx)*exps0);
    }
  }
  Width *= alphaP*pMass;

  //JMQ 190709 fix on Rb and  geometrical cross sections according to 
  //           Furihata's paper (JAERI-Data/Code 2001-105, p6)
  G4double Rb = 0.0;
  if (theA > 4) {
    Rb = 1.12*(resA13 + A13) - 0.86*((resA13 + A13)/(resA13*A13))+2.85;
  } else if (theA > 1) {
    Rb=1.5*(resA13 + A13);
  } else {
    Rb = 1.5*resA13;
  }

  G4double ild;
  if (exc < ExCN ) {
    G4double E0CN = ExCN - TCN*(G4Log(TCN) - 0.25*G4Log(a0) 
				- 1.25*G4Log(UxCN) 
				+ 2.0*std::sqrt(a0*UxCN));
    ild = G4Exp((exc-E0CN)/TCN)/TCN;
  } else {
    G4double x  = exc - delta0;
    G4double x1 = std::sqrt(a0*x);
    ild = G4Exp(2*x1)/(x*std::sqrt(x1));
  }

  Width *= (Rb*Rb/ild); 
  return Width;
}

G4Fragment* G4GEMProbabilityVI::Sample2DDistribution()
{
  G4Fragment* aFragment = nullptr;
  return aFragment;
}

G4double G4GEMProbabilityVI::I0(G4double t)
{
  return G4Exp(t) - 1.0;
}

G4double G4GEMProbabilityVI::I1(G4double t, G4double tx)
{
  return (t - tx + 1.0)*G4Exp(tx) - t - 1.0;
}

G4double G4GEMProbabilityVI::I2(G4double s0, G4double sx)
{
  G4double S = 1.0/std::sqrt(s0);
  G4double Sx = 1.0/std::sqrt(sx);
  
  G4double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
  G4double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*G4Exp(sx-s0);
  
  return p1-p2;
}

G4double G4GEMProbabilityVI::I3(G4double s0, G4double sx)
{
  G4double s2 = s0*s0;
  G4double sx2 = sx*sx;
  G4double S = 1.0/std::sqrt(s0);
  G4double S2 = S*S;
  G4double Sx = 1.0/std::sqrt(sx);
  G4double Sx2 = Sx*Sx;
  
  G4double p1 = S *(2.0 + S2 *( 4.0 + S2 *( 13.5 + S2 *( 60.0 + S2 * 325.125 ))));
  G4double p2 = Sx*Sx2 *((s2-sx2) + Sx2 *((1.5*s2+0.5*sx2) 
	      + Sx2 *((3.75*s2+0.25*sx2) + Sx2 *((12.875*s2+0.625*sx2) 
	      + Sx2 *((59.0625*s2+0.9375*sx2) + Sx2 *(324.8*s2+3.28*sx2))))));
  p2 *= G4Exp(sx-s0);
  return p1-p2; 
}

