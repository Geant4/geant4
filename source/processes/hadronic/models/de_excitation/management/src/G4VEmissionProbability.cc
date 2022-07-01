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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4VEmissionProbability::G4VEmissionProbability(G4int Z, G4int A)
  : OPTxs(3), pVerbose(1), theZ(Z), theA(A), elimit(CLHEP::MeV)
{
  pNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  pG4pow = G4Pow::GetInstance();
  if(A > 0) { pEvapMass = G4NucleiProperties::GetNuclearMass(theA, theZ); }
}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = pNuclearLevelData->GetParameters();
  OPTxs = param->GetDeexModelType();
  pVerbose = param->GetVerbose();
  fFD = param->GetDiscreteExcitationFlag();
  pTolerance = param->GetMinExcitation();
}

void G4VEmissionProbability::ResetIntegrator(size_t, G4double de, G4double eps)
{
  if(de > 0.0)  { elimit = de; }
  if(eps > 0.0) { accuracy = eps; }
}

G4double G4VEmissionProbability::EmissionProbability(const G4Fragment&, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::ComputeProbability(G4double, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::IntegrateProbability(G4double elow, 
                                                      G4double ehigh, 
                                                      G4double cb)
{
  pProbability = 0.0;
  if(elow >= ehigh) { return pProbability; }

  emin = elow;
  emax = ehigh;
  eCoulomb = cb;

  const G4double edeltamin = 0.2*CLHEP::MeV;
  const G4double edeltamax = 2*CLHEP::MeV;
  G4double edelta = std::max(std::min(elimit, edeltamax), edeltamin);
  G4double xbin = (emax - emin)/edelta + 1.0;
  G4int ibin = xbin;
  if(ibin < 4) ibin = 4;

  // providing smart binning 
  G4int nbin = ibin*5;
  edelta = (emax - emin)/ibin;

  G4double x(emin), y(0.0);
  G4double edelmicro = edelta*0.02;
  probmax = ComputeProbability(x + edelmicro, eCoulomb);
  G4double problast = probmax;
  if(pVerbose > 1) {
    G4cout << "### G4VEmissionProbability::IntegrateProbability: "
	   << "probmax=" << probmax << " Emin=" << emin
	   << " Emax=" << emax << " QB=" << cb << " nbin=" << nbin 
	   << G4endl;
  }
  fE1 = fE2 = fP2 = 0.0;
  G4double emax0 = emax - edelmicro; 
  G4bool endpoint = false;
  for(G4int i=0; i<nbin; ++i) {
    x += edelta;
    if(x >= emax0) { 
      x = emax0;
      endpoint = true;
    }
    y = ComputeProbability(x, eCoulomb);
    if(pVerbose > 2) { 
      G4cout << "    " << i << ".  E= " << x << "  prob= " << y
	     << " Edel= " << edelta << G4endl;
    } 
    if(y >= probmax) {
      probmax = y;
    } else if(0.0 == fE1 && 2*y < probmax) {
      fE1 = x;
    }
    
    G4double del = (y + problast)*edelta*0.5;
    pProbability += del;
    // end of the loop
    if(del < accuracy*pProbability || endpoint) { break; }
    problast = y;

    // smart step definition
    if(del != pProbability && del > 0.8*pProbability && 
       0.7*edelta > edeltamin) { 
      edelta *= 0.7;
    } else if(del < 0.1*pProbability && 1.5*edelta < edeltamax) { 
      edelta *= 1.5;
    }
  }
  if(fE1 > emin && fE1 < emax) {
    fE2 = std::max(0.5*(fE1 + emax), emax - edelta);
    fP2 = 2*ComputeProbability(fE2, eCoulomb);
  }

  if(pVerbose > 1) { 
    G4cout << " Probability= " << pProbability << " probmax= " 
           << probmax << " emin=" << emin << " emax=" << emax 
	   << " E1=" << fE1 << " E2=" << fE2 << G4endl; 
  }
  return pProbability;
}

G4double G4VEmissionProbability::SampleEnergy()
{
  static const G4double fact = 1.05;
  static const G4double alim = 0.05;
  static const G4double blim = 20.;
  probmax *= fact;

  // two regions with flat and exponential majorant 
  G4double del = emax - emin;
  G4double p1 = 1.0;
  G4double p2 = 0.0;
  G4double a0 = 0.0;
  G4double a1 = 1.0;
  G4double x;
  if(fE1 > 0.0 && fP2 > 0.0 && fP2 < 0.5*probmax) {
    a0 = G4Log(probmax/fP2)/(fE2 - fE1);
    del= fE1 - emin;
    p1 = del;
    x = a0*(emax - fE1);
    if(x < blim) {
      a1 = (x > alim) ? 1.0 - G4Exp(-x) : x*(1.0 - 0.5*x);
    }
    p2 = a1/a0;
    p1 /= (p1 + p2);
    p2 = 1.0 - p1;
  }  

  if(pVerbose > 1) {
    G4cout << "### G4VEmissionProbability::SampleEnergy: " 
	   << " Emin= " << emin << " Emax= " << emax 
           << "/n    E1=" << fE1 << " p1=" << p1 
	   << " probmax=" << probmax << " P2=" << fP2 << G4endl;
  }

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  const G4int nmax = 1000;
  G4double ekin, g, gmax;
  G4int n = 0;
  do {
    ++n;
    G4double q = rndm->flat();
    if(q <= p1) {
      gmax = probmax;
      ekin = del*q/p1 + emin;
    } else {
      ekin = fE1 - G4Log(1.0 - (q - p1)*a1/p2)/a0;
      x = a0*(ekin - fE1);
      gmax = fP2;
      if(x < blim) {
	gmax = probmax*((x > alim) ? G4Exp(-x) : 1.0 - x*(1.0 - 0.5*x));
      }
    }
    g = ComputeProbability(ekin, eCoulomb);
    if(pVerbose > 2) {
      G4cout << "    " << n
	     << ". prob= " << g << " probmax= " << probmax
	     << " Ekin= " << ekin << G4endl;
    }
    if((g > gmax || n > nmax) && pVerbose > 1) {
      G4cout << "### G4VEmissionProbability::SampleEnergy for Z= " << theZ 
             << " A= " << theA << " Eex(MeV)=" << fExc << " p1=" << p1
             << "\n    Warning n= " << n
	     << " prob/gmax=" << g/gmax 
	     << " prob=" << g << " gmax=" << gmax << " probmax=" << probmax 
	     << "\n    Ekin= " << ekin << " Emin= " << emin
	     << " Emax= " << emax << G4endl;
    }
  } while(gmax*rndm->flat() > g && n < nmax);
  G4double enew = FindRecoilExcitation(ekin);
  if(pVerbose > 1) {
    G4cout << "### SampleEnergy: Efinal= " 
	   << enew << " E=" << ekin << "  Eexc=" << fExcRes << G4endl;
  }
  return enew;
}

G4double G4VEmissionProbability::FindRecoilExcitation(const G4double e)
{
  G4double mass = pEvapMass + fExc;
    
  G4double m02 = pMass*pMass;
  G4double m12 = mass*mass;
  G4double m22 = pResMass*pResMass;
  G4double mres = std::sqrt(m02 + m12 - 2.*pMass*(mass + e));

  fExcRes = mres - pResMass;

  if(pVerbose > 1) {
    G4cout << "### FindRecoilExcitation for resZ= " 
           << resZ << " resA= " << resA 
           << " evaporated Z= " << theZ << " A= " << theA
	   << " Ekin= " << e << " Eexc= " << fExcRes << G4endl;
  }

  // residual nucleus is in the ground state
  if(fExcRes < pTolerance) {
    fExcRes = 0.0;
    return std::max(0.5*(m02 + m12 - m22)/pMass - mass, 0.0);
  }
  if(!fFD) { return e; }
 
  // select final state excitation
  auto lManager = pNuclearLevelData->GetLevelManager(resZ, resA);
  if(nullptr == lManager) { return e; }

  // levels are not known
  if(fExcRes > lManager->MaxLevelEnergy() + pTolerance) { return e; }

  // find level
  G4double elevel = lManager->NearestLevelEnergy(fExcRes);

  // excited level
  if(pMass > mass + pResMass + elevel && 
     std::abs(elevel - fExcRes) <= pTolerance) {
    G4double massR = pResMass + elevel;
    G4double mr2 = massR*massR;
    fExcRes = elevel;
    return std::max(0.5*(m02 + m12 - mr2)/pMass - mass, 0.0);
  }
  return e;
}
