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
#include "G4VSIntegration.hh"

G4PreCompoundFragment::G4PreCompoundFragment(const G4ParticleDefinition* p,
					     G4VCoulombBarrier* aCoulBarrier)
  : G4VPreCompoundFragment(p, aCoulBarrier)
{}

G4double G4PreCompoundFragment::CrossSection(G4double ekin)
{
  // compute power once
  if (OPTxs > 1 && 0 < index && theResA != lastA) {
    lastA = theResA;
    muu = G4KalbachCrossSection::ComputePowerParameter(lastA, index);
  }
  if (OPTxs == 0) { 
    recentXS = GetOpt0(ekin);
  } else if (OPTxs == 1) {
    G4int Z = std::min(theResZ, ZMAXNUCLEARDATA);
    const G4double lim = 2*CLHEP::MeV;
    const G4double Kmin = 20*CLHEP::keV;
    G4double e = lowEnergyLimitMeV[Z];
    if (e == 0.0) { e = lim; }
    G4double K = std::max(ekin, Kmin);
    e = std::max(e, K);
    recentXS = fXSection->GetElementCrossSection(e, Z)/CLHEP::millibarn;
    recentXS *= GetAlpha()*(1.0 + GetBeta()/K);

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
