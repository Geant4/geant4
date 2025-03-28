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
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modified:
// 06.09.2008 JMQ Also external choice has been added for:
//               - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko cleanup
//

#include "G4PreCompoundFragment.hh"
#include "G4KalbachCrossSection.hh"
#include "G4ChatterjeeCrossSection.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4InterfaceToXS.hh"
#include "G4IsotopeList.hh"
#include "Randomize.hh"

G4PreCompoundFragment::G4PreCompoundFragment(const G4ParticleDefinition* p,
					     G4VCoulombBarrier* aCoulBarrier)
  : G4VPreCompoundFragment(p, aCoulBarrier)
{}

G4double G4PreCompoundFragment::CalcEmissionProbability(const G4Fragment& fr)
{
  theEmissionProbability = (Initialize(fr)) ?
    IntegrateEmissionProbability(theMinKinEnergy, theMaxKinEnergy, fr) : 0.0;
  /*  
  G4cout << "## G4PreCompoundFragment::CalcEmisProb "
         << "Zf= " << fr.GetZ_asInt()
	 << " Af= " << fr.GetA_asInt()
	 << " Elow= " << theMinKinEnergy
	 << " Eup= " << theMaxKinEnergy
	 << " prob= " << theEmissionProbability
	 << " index=" << index << " Z=" << theZ << " A=" << theA
	 << G4endl;
  */
  return theEmissionProbability;
}

G4double 
G4PreCompoundFragment::IntegrateEmissionProbability(G4double low, G4double up,
                                                    const G4Fragment& fr)
{  
  static const G4double den = 1.0/CLHEP::MeV;
  G4double del = (up - low);
  G4int nbins = del*den;
  nbins = std::max(nbins, 4);
  del /= static_cast<G4double>(nbins);
  G4double e = low + 0.5*del;
  probmax = ProbabilityDistributionFunction(e, fr);
  //G4cout << "    0. e= " << e << "  y= " << probmax << G4endl;

  G4double sum = probmax;
  for (G4int i=1; i<nbins; ++i) {
    e += del;

    G4double y = ProbabilityDistributionFunction(e, fr); 
    probmax = std::max(probmax, y);
    sum += y;
    if(y < sum*0.01) { break; }
    //G4cout <<"   "<<i<<". e= "<<e<<"  y= "<<y<<" sum= "<<sum<< G4endl;
  }
  sum *= del;
  //G4cout << "Evap prob: " << sum << " probmax= " << probmax << G4endl;
  return sum;
}

G4double G4PreCompoundFragment::CrossSection(G4double ekin)
{
  /*
  G4cout << "G4PreCompoundFragment::CrossSection OPTxs=" << OPTxs << " E=" << ekin
	 << " resZ=" << theResZ << " resA=" << theResA << " index=" << index
	 << " fXSection:" << fXSection << G4endl;
  */
  // compute power once
  if (OPTxs > 1 && 0 < index && theResA != lastA) {
    lastA = theResA;
    muu = G4KalbachCrossSection::ComputePowerParameter(lastA, index);
  }
  if (OPTxs == 0) { 
    recentXS = GetOpt0(ekin);
  } else if (OPTxs == 1) {
    G4int Z = std::min(theResZ, ZMAXNUCLEARDATA);
    //G4double e = std::max(ekin, lowEnergyLimitMeV[Z]);
    recentXS = fXSection->GetElementCrossSection(ekin, Z)/CLHEP::millibarn;

  } else if (OPTxs == 2) { 
    recentXS = G4ChatterjeeCrossSection::ComputeCrossSection(ekin, 
                                                             theCoulombBarrier, 
							     theResA13, muu, 
							     index, theZ, theResA); 

  } else { 
    recentXS = G4KalbachCrossSection::ComputeCrossSection(ekin, theCoulombBarrier, 
						          theResA13, muu, index,
						          theZ, theA, theResA);
  }
  return recentXS;
}  

G4double G4PreCompoundFragment::GetOpt0(G4double ekin) const
// OPT=0 : Dostrovski's cross section
{
  G4double r0 = theParameters->GetR0()*theResA13;
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  // with the rest of the options
  return 1.e+25*CLHEP::pi*r0*r0*theResA13*GetAlpha()*(1.0 + GetBeta()/ekin);
}

G4double G4PreCompoundFragment::SampleKineticEnergy(const G4Fragment& fragment) 
{
  G4double delta = theMaxKinEnergy - theMinKinEnergy;
  static const G4double toler = 1.25;
  probmax *= toler;
  G4double prob, T(0.0);
  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  G4int i;
  for(i=0; i<100; ++i) {
    T = theMinKinEnergy + delta*rndm->flat();
    prob = ProbabilityDistributionFunction(T, fragment);
    /*
    if(prob > probmax) { 
      G4cout << "G4PreCompoundFragment WARNING: prob= " << prob 
	     << " probmax= " << probmax << G4endl;
      G4cout << "i= " << i << " Z= " << theZ << " A= " << theA 
	     << " resZ= " << theResZ << " resA= " << theResA << "\n"
	     << " T= " << T << " Tmax= " << theMaxKinEnergy 
	     << " Tmin= " << limit
	     << G4endl;
      for(G4int i=0; i<N; ++i) { G4cout << " " << probability[i]; }
      G4cout << G4endl; 
    }
    */
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    if(probmax*rndm->flat() <= prob) { break; }
  }
  /*
  G4cout << "G4PreCompoundFragment: i= " << i << " Z= " << theZ 
         << " A= " << theA <<"  T(MeV)= " << T << " Emin(MeV)= " 
         << theMinKinEnergy << " Emax= " << theMaxKinEnergy << G4endl;
  */
  return T;
}

