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
// $Id: G4VPreCompoundFragment.cc,v 1.13 2010-08-28 15:16:55 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// J. M. Quesada (August 2008).  Based  on previous work by V. Lara
//
// Modified:
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup

#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4NucleiProperties.hh"

G4VPreCompoundFragment::G4VPreCompoundFragment(
  const G4ParticleDefinition* part, G4VCoulombBarrier* aCoulombBarrier)
  : particle(part), theCoulombBarrierPtr(aCoulombBarrier),
    theRestNucleusA(0),theRestNucleusZ(0),theBindingEnergy(0.0), 
    theMaximalKineticEnergy(-MeV),theRestNucleusMass(0.0),
    theReducedMass(0.0),theMomentum(0.,0.,0.,0.),
    theEmissionProbability(0.0),theCoulombBarrier(0.0)
{
  theA = particle->GetBaryonNumber();
  theZ = G4int(particle->GetPDGCharge()/eplus + 0.1);
  theMass = particle->GetPDGMass();
  theParameters = G4PreCompoundParameters::GetAddress();
  g4pow = G4Pow::GetInstance();
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
  std::ios::fmtflags old_floatfield = out.flags();
  out.setf(std::ios::floatfield);
    
  out 
    << "PreCompoundModel Emitted Fragment: A = " 
    << std::setprecision(3) << theFragment->theA 
    << ", Z = " << std::setprecision(3) << theFragment->theZ;
    out.setf(std::ios::scientific,std::ios::floatfield);
    //   out
    //     << ", U = " << theFragment->theExcitationEnergy/MeV 
    //     << " MeV" << endl
    //     << "          P = (" 
    //     << theFragment->theMomentum.x()/MeV << ","
    //     << theFragment->theMomentum.y()/MeV << ","
    //     << theFragment->theMomentum.z()/MeV 
    //     << ") MeV   E = " 
    //     << theFragment->theMomentum.t()/MeV << " MeV";
    
    out.setf(old_floatfield,std::ios::floatfield);
    return out;
}

void G4VPreCompoundFragment::
Initialize(const G4Fragment & aFragment)
{
  theRestNucleusA = aFragment.GetA_asInt() - theA;
  theRestNucleusZ = aFragment.GetZ_asInt() - theZ;
  theRestNucleusA13 = g4pow->Z13(theRestNucleusA);

  if ((theRestNucleusA < theRestNucleusZ) ||
      (theRestNucleusA < theA) ||
      (theRestNucleusZ < theZ)) 
    {
      // In order to be sure that emission probability will be 0.
      theMaximalKineticEnergy = 0.0;
      return;
    }
    
  // Calculate Coulomb barrier
  theCoulombBarrier = theCoulombBarrierPtr->
    GetCoulombBarrier(theRestNucleusA,theRestNucleusZ,
		      aFragment.GetExcitationEnergy());

  // Calculate masses
  theRestNucleusMass = 
    G4NucleiProperties::GetNuclearMass(theRestNucleusA, theRestNucleusZ);
  theReducedMass = theRestNucleusMass*theMass/(theRestNucleusMass + theMass);

  // Compute Binding Energies for fragments 
  // (needed to separate a fragment from the nucleus)
  theBindingEnergy = theRestNucleusMass + theMass - aFragment.GetGroundStateMass();
  
  //theBindingEnergy = G4NucleiProperties::GetMassExcess(static_cast<G4int>(theA),static_cast<G4int>(theZ)) +
  //G4NucleiProperties::GetMassExcess(static_cast<G4int>(theRestNucleusA),static_cast<G4int>(theRestNucleusZ)) -
  //G4NucleiProperties::GetMassExcess(static_cast<G4int>(aFragment.GetA()),static_cast<G4int>(aFragment.GetZ()));
  
  // Compute Maximal Kinetic Energy which can be carried by fragments after separation
  // This is the true (assimptotic) maximal kinetic energy
  G4double m  = aFragment.GetMomentum().m();
  G4double rm = theRestNucleusMass;
  G4double em = theMass;
  theMaximalKineticEnergy = ((m - rm)*(m + rm) + em*em)/(2.0*m) - em;
 
  return;
}
