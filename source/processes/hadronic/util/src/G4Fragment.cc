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
// 27.10.2021 A.Ribon extension for hypernuclei.

#include "G4Fragment.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4ios.hh"
#include <iomanip>

G4Allocator<G4Fragment>*& pFragmentAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4Fragment>* _instance = nullptr;
  return _instance;
}

// Default constructor
G4Fragment::G4Fragment() :
  theA(0),
  theZ(0),
  theL(0),
  theExcitationEnergy(0.0),
  theGroundStateMass(0.0),
  theMomentum(G4LorentzVector(0,0,0,0)),
  thePolarization(nullptr),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  xLevel(0),
  theParticleDefinition(nullptr),
  spin(0.0),
  theCreationTime(0.0)
{}

// Copy Constructor
G4Fragment::G4Fragment(const G4Fragment &right) :
   theA(right.theA),
   theZ(right.theZ),
   theL(right.theL),
   theExcitationEnergy(right.theExcitationEnergy),
   theGroundStateMass(right.theGroundStateMass),
   theMomentum(right.theMomentum),
   thePolarization(right.thePolarization),
   creatorModel(right.creatorModel),
   numberOfParticles(right.numberOfParticles),
   numberOfCharged(right.numberOfCharged),
   numberOfHoles(right.numberOfHoles),
   numberOfChargedHoles(right.numberOfChargedHoles),
   numberOfShellElectrons(right.numberOfShellElectrons),
   xLevel(right.xLevel),
   theParticleDefinition(right.theParticleDefinition),
   spin(right.spin),
   theCreationTime(right.theCreationTime),
   isLongLived(right.isLongLived)
{}

G4Fragment::~G4Fragment()
{}

G4Fragment::G4Fragment(G4int A, G4int Z, const G4LorentzVector& aMomentum) :
  theA(A),
  theZ(Z),
  theL(0),
  theExcitationEnergy(0.0),
  theGroundStateMass(0.0),
  theMomentum(aMomentum),
  thePolarization(nullptr),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  xLevel(0),
  theParticleDefinition(nullptr),
  spin(0.0),
  theCreationTime(0.0)
{
  if(theA > 0) { 
    CalculateMassAndExcitationEnergy();
  }
}

G4Fragment::G4Fragment(G4int A, G4int Z, G4int numberOfLambdas,
                       const G4LorentzVector& aMomentum) :
  theA(A),
  theZ(Z),
  theL(std::max(numberOfLambdas,0)),
  theExcitationEnergy(0.0),
  theGroundStateMass(0.0),
  theMomentum(aMomentum),
  thePolarization(nullptr),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  xLevel(0),
  theParticleDefinition(nullptr),
  spin(0.0),
  theCreationTime(0.0)
{
  if(theA > 0) { 
    CalculateMassAndExcitationEnergy();
  }
}

// This constructor is for initialize photons or electrons
G4Fragment::G4Fragment(const G4LorentzVector& aMomentum, 
		       const G4ParticleDefinition * aParticleDefinition) :
  theA(0),
  theZ(0),
  theL(0),
  theExcitationEnergy(0.0),
  theMomentum(aMomentum),
  thePolarization(nullptr),
  creatorModel(-1),
  numberOfParticles(0),
  numberOfCharged(0),
  numberOfHoles(0),
  numberOfChargedHoles(0),
  numberOfShellElectrons(0),
  xLevel(0),
  theParticleDefinition(aParticleDefinition),
  spin(0.0),
  theCreationTime(0.0)
{
  if(aParticleDefinition->GetPDGEncoding() != 22 && 
     aParticleDefinition->GetPDGEncoding() != 11) {
    G4String text = "G4Fragment::G4Fragment constructor for gamma used for "
      + aParticleDefinition->GetParticleName();  
    throw G4HadronicException(__FILE__, __LINE__, text);
  }
  theGroundStateMass = aParticleDefinition->GetPDGMass();
}

void G4Fragment::CalculateMassAndExcitationEnergy()
{
  // check input
  if(theZ > theA || theZ + theL > theA) {
    G4String text = "G4Fragment::CalculateMassAndExcitationEnergy: inconsistent number of nucleons is ignored";
    G4cout << text << G4endl; 
    G4cout << "       Z=" << theZ << " A=" << theA 
	   << " nLambdas=" << theL << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, text);
  }
  // compute mass
  theGroundStateMass = ( theL == 0 )
    ? G4NucleiProperties::GetNuclearMass(theA, theZ)
    : G4HyperNucleiProperties::GetNuclearMass(theA, theZ, theL); 

  // excitation energy
  const G4double minFragExcitation = 10.*CLHEP::eV;
  theExcitationEnergy = theMomentum.mag() - theGroundStateMass;
  if(theExcitationEnergy < minFragExcitation) {
    if(theExcitationEnergy < -minFragExcitation) {
      ExcitationEnergyWarning();
    }
    theExcitationEnergy = 0.0;
  }
}

void G4Fragment::SetExcEnergyAndMomentum(G4double eexc, 
					 const G4LorentzVector& v)
{
  theExcitationEnergy = eexc;
  theMomentum.set(0.0, 0.0, 0.0, theGroundStateMass + eexc);
  theMomentum.boost(v.boostVector());
}

G4double G4Fragment::GetBindingEnergy() const
{
  const G4double lambdaMass = 1.115683*CLHEP::GeV;
  return (theA-theZ-theL)*CLHEP::neutron_mass_c2 
    + theZ*CLHEP::proton_mass_c2 + theL*lambdaMass
    - theGroundStateMass;
}

G4Fragment & G4Fragment::operator=(const G4Fragment &right)
{
  if (this != &right) {
    theA = right.theA;
    theZ = right.theZ;
    theL = right.theL;
    theExcitationEnergy = right.theExcitationEnergy;
    theGroundStateMass = right.theGroundStateMass;
    theMomentum  = right.theMomentum;
    thePolarization = right.thePolarization;
    creatorModel = right.creatorModel;
    numberOfParticles = right.numberOfParticles;
    numberOfCharged = right.numberOfCharged;
    numberOfHoles = right.numberOfHoles;
    numberOfChargedHoles = right.numberOfChargedHoles;
    numberOfShellElectrons = right.numberOfShellElectrons;
    xLevel = right.xLevel;
    theParticleDefinition = right.theParticleDefinition;
    spin = right.spin;
    theCreationTime = right.theCreationTime;
    isLongLived = right.isLongLived;
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

std::ostream& operator << (std::ostream &out, const G4Fragment &theFragment)
{
  std::ios::fmtflags old_floatfield = out.flags();
  out.setf(std::ios::floatfield);

  out << "Fragment: A = " << std::setw(3) << theFragment.theA 
      << ", Z = " << std::setw(3) << theFragment.theZ
      << ", numberOfLambdas = " << std::setw(3) << theFragment.theL ;
  out.setf(std::ios::scientific,std::ios::floatfield);

  // Store user's precision setting and reset to (3) here: back-compatibility
  std::streamsize floatPrec = out.precision();

  out << std::setprecision(3)
      << ", U = " << theFragment.GetExcitationEnergy()/CLHEP::MeV 
      << " MeV  ";
  if(theFragment.GetCreatorModelID() >= 0) { 
    out << " creatorModelID= " << theFragment.GetCreatorModelID(); 
  }
  if(theFragment.GetCreationTime() > 0.0) { 
    out << "  Time= " << theFragment.GetCreationTime()/CLHEP::ns << " ns"; 
  }
  out << G4endl
      << "          P = (" 
      << theFragment.GetMomentum().x()/CLHEP::MeV << ","
      << theFragment.GetMomentum().y()/CLHEP::MeV << ","
      << theFragment.GetMomentum().z()/CLHEP::MeV 
      << ") MeV   E = " 
      << theFragment.GetMomentum().t()/CLHEP::MeV << " MeV"
      << G4endl;

  out << "    #spin= " << theFragment.GetSpin()
      << "    #floatLevelNo= " << theFragment.GetFloatingLevelNumber() << "  ";
        
  if (theFragment.GetNumberOfExcitons() != 0) {
    out << "   " 
	<< "#Particles= " << theFragment.GetNumberOfParticles() 
	<< ", #Charged= " << theFragment.GetNumberOfCharged()
	<< ", #Holes= "   << theFragment.GetNumberOfHoles()
	<< ", #ChargedHoles= " << theFragment.GetNumberOfChargedHoles();
  } 
  out << G4endl;
  if(theFragment.GetNuclearPolarization()) { 
    out << *(theFragment.GetNuclearPolarization()); 
  }
  //out << G4endl;
  out.setf(old_floatfield,std::ios::floatfield);
  out.precision(floatPrec);

  return out;
}

void G4Fragment::ExcitationEnergyWarning()
{
#ifdef G4VERBOSE
  G4cout << "G4Fragment::CalculateExcitationEnergy(): WARNING "
	 << " GraundStateMass(MeV)= " << theGroundStateMass 
	 <<G4endl;
  G4cout << *this << G4endl;
#endif
}

void G4Fragment::NumberOfExitationWarning(const G4String& value)
{
  G4cout << "G4Fragment::"<< value << " ERROR "
	 << G4endl;
  G4cout << this << G4endl; 
  G4String text = "G4Fragment::G4Fragment wrong exciton number ";
  throw G4HadronicException(__FILE__, __LINE__, text);
}

void G4Fragment::SetAngularMomentum(const G4ThreeVector& v)
{
  spin = v.mag();
}

G4ThreeVector G4Fragment::GetAngularMomentum() const
{
  G4ThreeVector v(0.0,0.0,spin);
  return v;
}
