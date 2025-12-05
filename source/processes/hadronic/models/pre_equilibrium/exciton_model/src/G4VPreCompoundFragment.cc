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
// J. M. Quesada (August 2008).  Based  on previous work by V. Lara
//
// Modified:
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup

#include "G4VPreCompoundFragment.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4VCoulombBarrier.hh"
#include "G4InterfaceToXS.hh"

G4VPreCompoundFragment::G4VPreCompoundFragment(
  const G4ParticleDefinition* part, G4VCoulombBarrier* aCoulombBarrier)
  : theA(part->GetBaryonNumber()),
    theZ(G4lrint(part->GetPDGCharge()/CLHEP::eplus)),
    particle(part),
    theCoulombBarrierPtr(aCoulombBarrier)
{
  theMass = particle->GetPDGMass();
  fNucData = G4NuclearLevelData::GetInstance();
  theParameters = fNucData->GetParameters();
  OPTxs = theParameters->GetDeexModelType();
  g4calc = G4Pow::GetInstance();

  if (1 == theZ && 1 == theA) { index = 1; }
  else if (1 == theZ && 2 == theA) { index = 2; }
  else if (1 == theZ && 3 == theA) { index = 3; }
  else if (2 == theZ && 3 == theA) { index = 4; }
  else if (2 == theZ && 4 == theA) { index = 5; }

  if (OPTxs == 1) {
    fXSection = new G4InterfaceToXS(particle, index);
  }
  InitialiseIntegrator(0.005, 0.25, 1.10, 0.5*CLHEP::MeV,
		       0.2*CLHEP::MeV, 5*CLHEP::MeV);
}

G4VPreCompoundFragment::~G4VPreCompoundFragment()
{
  delete theCoulombBarrierPtr;
  delete fXSection;
}

std::ostream& 
operator << (std::ostream &out, const G4VPreCompoundFragment &theFragment)
{
  out << &theFragment;
  return out; 
}

std::ostream& 
operator << (std::ostream &out, const G4VPreCompoundFragment *theFragment)
{
  out 
    << "PreCompoundModel Emitted Fragment: Z= " << theFragment->GetZ() 
    << " A= " << theFragment->GetA()
    << " Mass(GeV)= " << theFragment->GetNuclearMass()/CLHEP::GeV;
  return out;
}

G4bool 
G4VPreCompoundFragment::Initialize(const G4Fragment& aFragment)
{
  theFragA = aFragment.GetA_asInt();
  theFragZ = aFragment.GetZ_asInt();
  theResA = theFragA - theA;
  theResZ = theFragZ - theZ;

  theMinKinEnergy = theMaxKinEnergy = theCoulombBarrier = 0.0;
  if ((theResA < theResZ) || (theResA < theA) || (theResZ < theZ)
      || (theResA == theA && theResZ < theZ)
      || ((theResA > 1) && (theResA == theResZ || theResZ == 0))) {
    return false;
  }
  pFragment = &aFragment;
  theResMass = G4NucleiProperties::GetNuclearMass(theResA, theResZ);
  G4double Ecm = aFragment.GetMomentum().m();
  if (Ecm <= theResMass + theMass) { return 0.0; }

  theResA13 = g4calc->Z13(theResA);
  G4int nex = aFragment.GetNumberOfExcitons();

  G4double elim = 0.0;
  if (0 < theZ) {
    theCoulombBarrier = theCoulombBarrierPtr->
      GetCoulombBarrier(theResA + nex, theResZ, aFragment.GetExcitationEnergy());
    elim = (0 < OPTxs) ? theCoulombBarrier*0.5 : theCoulombBarrier;
  }
      
  // Compute Maximal Kinetic Energy which can be carried by fragments 
  // after separation - the true assimptotic value
  theMaxKinEnergy =
    0.5*((Ecm - theResMass)*(Ecm + theResMass) + theMass*theMass)/Ecm - theMass;
  G4double resM = Ecm - theMass - elim;
  if (resM < theResMass) { return false; }
  theMinKinEnergy =
    0.5*((Ecm - resM)*(Ecm + resM) + theMass*theMass)/Ecm - theMass;
  theMinKinEnergy = std::max(theMinKinEnergy, 0.0);

  if (theMinKinEnergy >= theMaxKinEnergy) { return false; }
  // Calculate masses
  theReducedMass = theResMass*theMass/(theResMass + theMass);

  // Compute Binding Energies for fragments 
  // needed to separate a fragment from the nucleus
  theBindingEnergy = theResMass + theMass - aFragment.GetGroundStateMass();
  return true;
}

G4double G4VPreCompoundFragment::CalcEmissionProbability(const G4Fragment& fr)
{
  G4bool ok = Initialize(fr);
  theEmissionProbability = 0.0;
  if (ok) {
    theEmissionProbability = ComputeIntegral(theMinKinEnergy, theMaxKinEnergy);
  }
  return theEmissionProbability;
}

G4double G4VPreCompoundFragment::SampleKineticEnergy(const G4Fragment&)
{
  G4double ekin = SampleValue();
  return ekin;
}

G4double G4VPreCompoundFragment::ProbabilityDensityFunction(G4double ekin)
{
  G4double e = std::max(ekin, 0.02);
  return ProbabilityDistributionFunction(e, *pFragment);
}
