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
// $Id: G4Fragment.cc 85824 2014-11-05 15:26:17Z gcosmo $
//
//---------------------------------------------------------------------
//
// Geant4 class G4Fragment
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
// Modifications:
// 03.05.2010 V.Ivanchenko General cleanup; moved obsolete methods from
//            inline to source 
// 25.09.2010 M. Kelsey -- Change "setprecision" to "setwidth" in printout,
//	      add null pointer check.

#include "G4Fragment.hh"
#include "G4HadronicException.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ios.hh"
#include <iomanip>

//#define debug_G4Fragment 

G4ThreadLocal G4Allocator<G4Fragment> *pFragmentAllocator = 0;

// Default constructor
G4Fragment::G4Fragment() :
  theA(0),
  theZ(0),
  theExcitationEnergy(0.0),
  theGroundStateMass(0.0),
  theMomentum(G4LorentzVector(0,0,0,0)),
  theAngularMomentum(G4ThreeVector(0,0,0)),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  theParticleDefinition(0),
  theCreationTime(0.0),
  isStable(true)
{}

// Copy Constructor
G4Fragment::G4Fragment(const G4Fragment &right) 
{
   theA = right.theA;
   theZ = right.theZ;
   theExcitationEnergy = right.theExcitationEnergy;
   theGroundStateMass = right.theGroundStateMass;
   theMomentum  = right.theMomentum;
   theAngularMomentum = right.theAngularMomentum;
   creatorModel = right.creatorModel;
   numberOfParticles = right.numberOfParticles;
   numberOfCharged = right.numberOfCharged;
   numberOfHoles = right.numberOfHoles;
   numberOfChargedHoles = right.numberOfChargedHoles;
   numberOfShellElectrons = right.numberOfShellElectrons;
   theParticleDefinition = right.theParticleDefinition;
   theCreationTime = right.theCreationTime;
   isStable = right.isStable;
}

G4Fragment::~G4Fragment()
{}

G4Fragment::G4Fragment(G4int A, G4int Z, const G4LorentzVector& aMomentum) :
  theA(A),
  theZ(Z),
  theMomentum(aMomentum),
  theAngularMomentum(G4ThreeVector(0,0,0)),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  theParticleDefinition(0),
  theCreationTime(0.0),
  isStable(true)
{
  theExcitationEnergy = 0.0;
  theGroundStateMass = 0.0;
  if(theA > 0) { 
    CalculateGroundStateMass();
    CalculateExcitationEnergy(); 

    // default flag of stability for excited fragments is false
    // it may be overwritten by SetStable(G4bool val) method
    if(theExcitationEnergy > 0.0) { isStable = false; }
  }
}

// This constructor is for initialize photons or electrons
G4Fragment::G4Fragment(const G4LorentzVector& aMomentum, 
		       const G4ParticleDefinition * aParticleDefinition) :
  theA(0),
  theZ(0),
  theMomentum(aMomentum),
  theAngularMomentum(G4ThreeVector(0,0,0)),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  theParticleDefinition(aParticleDefinition),
  theCreationTime(0.0),
  isStable(true)
{
  theExcitationEnergy = 0.0;
  if(aParticleDefinition->GetPDGEncoding() != 22 && 
     aParticleDefinition->GetPDGEncoding() != 11) {
    G4String text = "G4Fragment::G4Fragment constructor for gamma used for "
      + aParticleDefinition->GetParticleName();  
    throw G4HadronicException(__FILE__, __LINE__, text);
  }
  theGroundStateMass = aParticleDefinition->GetPDGMass();
}

G4Fragment & G4Fragment::operator=(const G4Fragment &right)
{
  if (this != &right) {
    theA = right.theA;
    theZ = right.theZ;
    theExcitationEnergy = right.theExcitationEnergy;
    theGroundStateMass = right.theGroundStateMass;
    theMomentum  = right.theMomentum;
    theAngularMomentum = right.theAngularMomentum;
    creatorModel = right.creatorModel;
    numberOfParticles = right.numberOfParticles;
    numberOfCharged = right.numberOfCharged;
    numberOfHoles = right.numberOfHoles;
    numberOfChargedHoles = right.numberOfChargedHoles;
    numberOfShellElectrons = right.numberOfShellElectrons;
    theParticleDefinition = right.theParticleDefinition;
    theCreationTime = right.theCreationTime;
    isStable = right.isStable;
  }
  return *this;
}

G4bool G4Fragment::operator==(const G4Fragment &right) const
{
  return (this == (G4Fragment *) &right);
}

G4bool G4Fragment::operator!=(const G4Fragment &right) const
{
  return (this != (G4Fragment *) &right);
}

std::ostream& operator << (std::ostream &out, const G4Fragment *theFragment)
{
  if (!theFragment) {
    out << "Fragment: null pointer ";
    return out;
  }

  std::ios::fmtflags old_floatfield = out.flags();
  out.setf(std::ios::floatfield);

  out << "Fragment: A = " << std::setw(3) << theFragment->theA 
      << ", Z = " << std::setw(3) << theFragment->theZ ;
  out.setf(std::ios::scientific,std::ios::floatfield);

  // Store user's precision setting and reset to (3) here: back-compatibility
  std::streamsize floatPrec = out.precision();

  out << std::setprecision(3)
      << ", U = " << theFragment->GetExcitationEnergy()/CLHEP::MeV 
      << " MeV  IsStable= " << theFragment->IsStable();
  if(theFragment->creatorModel >= 0) { 
    out << " creatorModelType= " << theFragment->creatorModel; 
  }
  out << G4endl
      << "          P = (" 
      << theFragment->theMomentum.x()/CLHEP::MeV << ","
      << theFragment->theMomentum.y()/CLHEP::MeV << ","
      << theFragment->theMomentum.z()/CLHEP::MeV 
      << ") MeV   E = " 
      << theFragment->theMomentum.t()/CLHEP::MeV << " MeV"
      << G4endl;
       
  // What about Angular momentum???
  if (theFragment->GetNumberOfExcitons() != 0) {
    out << "          " 
	<< "#Particles= " << theFragment->numberOfParticles 
	<< ", #Charged= " << theFragment->numberOfCharged
	<< ", #Holes= "   << theFragment->numberOfHoles
	<< ", #ChargedHoles= " << theFragment->numberOfChargedHoles
	<< G4endl;
  }
  out.setf(old_floatfield,std::ios::floatfield);
  out.precision(floatPrec);

  return out;
}

std::ostream& operator << (std::ostream &out, const G4Fragment &theFragment)
{
  out << &theFragment;
  return out; 
}

void G4Fragment::ExcitationEnergyWarning()
{
  if (theExcitationEnergy < -10 * CLHEP::eV) {

#ifdef G4VERBOSE
    G4cout << "G4Fragment::CalculateExcitationEnergy(): WARNING "<<G4endl;
    G4cout << *this << G4endl;
#endif

#ifdef debug_G4Fragment
    G4ExceptionDescription ed;
    ed << *this << G4endl;
    G4Exception("G4Fragment::ExcitationEnergyWarning()", "had777", 
		FatalException,ed);
#endif
  }
  theExcitationEnergy = 0.0;
}

void G4Fragment::NumberOfExitationWarning(const G4String& value)
{
  G4cout << "G4Fragment::"<< value << " ERROR "
	 << G4endl;
  G4cout << this << G4endl; 
  G4String text = "G4Fragment::G4Fragment wrong exciton number ";
  throw G4HadronicException(__FILE__, __LINE__, text);
}
