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
// $Id$
//
// J. M. Quesada (August 2008).  Based  on previous work by V. Lara
//
// Modified:
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup

#include "G4VPreCompoundFragment.hh"
#include "G4SystemOfUnits.hh"
#include "G4PreCompoundParameters.hh"
#include "G4NucleiProperties.hh"

G4VPreCompoundFragment::G4VPreCompoundFragment(
  const G4ParticleDefinition* part, G4VCoulombBarrier* aCoulombBarrier)
  : particle(part), theCoulombBarrierPtr(aCoulombBarrier),
    theRestNucleusA(0),theRestNucleusZ(0),theBindingEnergy(0.0), 
    theMaximalKineticEnergy(-MeV),theRestNucleusMass(0.0),
    theReducedMass(0.0),theMomentum(0.,0.,0.,0.),
    theEmissionProbability(0.0),theCoulombBarrier(0.0),
    OPTxs(3),useSICB(false)
{
  theA = particle->GetBaryonNumber();
  theZ = G4int(particle->GetPDGCharge()/eplus + 0.1);
  theMass = particle->GetPDGMass();
  theParameters = G4PreCompoundParameters::GetAddress();
  g4pow = G4Pow::GetInstance();
  theRestNucleusA13 = 0;
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
  theRestNucleusA = aFragment.GetA_asInt() - theA;
  theRestNucleusZ = aFragment.GetZ_asInt() - theZ;

  if ((theRestNucleusA < theRestNucleusZ) ||
      (theRestNucleusA < theA) ||
      (theRestNucleusZ < theZ)) 
    {
      // In order to be sure that emission probability will be 0.
      theMaximalKineticEnergy = 0.0;
      return;
    }

  theRestNucleusA13 = g4pow->Z13(theRestNucleusA);
    
  // Calculate Coulomb barrier
  theCoulombBarrier = theCoulombBarrierPtr->
    GetCoulombBarrier(theRestNucleusA,theRestNucleusZ,
		      aFragment.GetExcitationEnergy());

  // Calculate masses
  theRestNucleusMass = 
    G4NucleiProperties::GetNuclearMass(theRestNucleusA, theRestNucleusZ);
  theReducedMass = theRestNucleusMass*theMass/(theRestNucleusMass + theMass);

  // Compute Binding Energies for fragments 
  // needed to separate a fragment from the nucleus
  theBindingEnergy = 
    theRestNucleusMass + theMass - aFragment.GetGroundStateMass();
    
  // Compute Maximal Kinetic Energy which can be carried by fragments 
  // after separation - the true assimptotic value
  G4double Ecm  = aFragment.GetMomentum().m();
  theMaximalKineticEnergy = 
    ((Ecm-theRestNucleusMass)*(Ecm+theRestNucleusMass) + theMass*theMass)/(2.0*Ecm)-theMass;
}
