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

G4VPreCompoundFragment::G4VPreCompoundFragment(
  const G4ParticleDefinition* part, G4VCoulombBarrier* aCoulombBarrier)
  : particle(part), theCoulombBarrierPtr(aCoulombBarrier),
    theMomentum(0.,0.,0.,0.),
    theA(particle->GetBaryonNumber()),
    theZ(G4lrint(particle->GetPDGCharge())),
    theResA(0),theResZ(0),theFragA(0),theFragZ(0),theBindingEnergy(0.0), 
    theMinKinEnergy(0.0),theMaxKinEnergy(0.0),theResMass(0.0),
    theReducedMass(0.0),
    theEmissionProbability(0.0),theCoulombBarrier(0.0),
    OPTxs(3),useSICB(true)
{
  theMass = particle->GetPDGMass();
  fNucData = G4NuclearLevelData::GetInstance();
  theParameters = fNucData->GetParameters();
  g4calc = G4Pow::GetInstance();
  theResA13 = 0.0;
}

G4VPreCompoundFragment::~G4VPreCompoundFragment()
{}

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

void 
G4VPreCompoundFragment::Initialize(const G4Fragment & aFragment)
{
  theFragA = aFragment.GetA_asInt();
  theFragZ = aFragment.GetZ_asInt();
  theResA = theFragA - theA;
  theResZ = theFragZ - theZ;

  theMinKinEnergy = theMaxKinEnergy = theCoulombBarrier = 0.0;
  if ((theResA < theResZ) || (theResA < theA) || (theResZ < theZ)) {
    return;
  }

  theResA13 = g4calc->Z13(theResA);
  theCoulombBarrier = theCoulombBarrierPtr->
    GetCoulombBarrier(theResA,theResZ,aFragment.GetExcitationEnergy());
    
  G4double elim = (0 == OPTxs) ? theCoulombBarrier : theCoulombBarrier*0.5;
  
  // Calculate masses
  theResMass = G4NucleiProperties::GetNuclearMass(theResA, theResZ);
  theReducedMass = theResMass*theMass/(theResMass + theMass);

  // Compute Binding Energies for fragments 
  // needed to separate a fragment from the nucleus
  theBindingEnergy = theResMass + theMass - aFragment.GetGroundStateMass();
    
  // Compute Maximal Kinetic Energy which can be carried by fragments 
  // after separation - the true assimptotic value
  G4double Ecm = aFragment.GetMomentum().m();
  G4double twoEcm = Ecm + Ecm;
  theMaxKinEnergy = std::max(((Ecm-theResMass)*(Ecm+theResMass) + theMass*theMass)
			     /twoEcm - theMass,0.0);
  theMinKinEnergy = (elim == 0.0) ? 0.0 :
    std::max(((theMass+elim)*(twoEcm-theMass-elim) + 
	      theMass*theMass)/twoEcm - theMass,0.0);
}
