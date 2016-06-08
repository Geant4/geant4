//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Fragment.cc,v 1.14 2002/06/10 13:27:54 jwellisc Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#include "G4Fragment.hh"


// Default constructor
G4Fragment::G4Fragment() :
  theA(0),
  theZ(0),
  theExcitationEnergy(0.0),
  theMomentum(0),
  theAngularMomentum(0),
  numberOfParticles(0),
  numberOfHoles(0),
  numberOfCharged(0),
  theParticleDefinition(0),
  theCreationTime(0.0)
#ifdef pctest
  ,theCreatorModel("No name")
#endif
{
  theAngularMomentum = IsotropicRandom3Vector();
}

// Copy Constructor
G4Fragment::G4Fragment(const G4Fragment &right) 
{
   theA = right.theA;
   theZ = right.theZ;
   theExcitationEnergy = right.theExcitationEnergy;
   theMomentum  = right.theMomentum;
   theAngularMomentum = right.theAngularMomentum;
   numberOfParticles = right.numberOfParticles;
   numberOfHoles = right.numberOfHoles;
   numberOfCharged = right.numberOfCharged;
   theParticleDefinition = right.theParticleDefinition;
   theCreationTime = right.theCreationTime;
#ifdef pctest
   theCreatorModel = right.theCreatorModel;
#endif
}


G4Fragment::~G4Fragment()
{
}


G4Fragment::G4Fragment(const G4int A, const G4int Z, const G4LorentzVector aMomentum) :
  theA(A),
  theZ(Z),
  theMomentum(aMomentum),
  numberOfParticles(0),
  numberOfHoles(0),
  numberOfCharged(0),
  theParticleDefinition(0),
  theCreationTime(0.0)
#ifdef pctest
  ,theCreatorModel("No name")
#endif
{
  theExcitationEnergy = theMomentum.mag() - 
                        G4ParticleTable::GetParticleTable()->GetIonTable()
			->GetIonMass( static_cast<G4int>(theZ), static_cast<G4int>(theA) );
  if( theExcitationEnergy < 0.0 )
    if( theExcitationEnergy > -10.0 * eV )
      theExcitationEnergy = 0.0;
    else 
    {
      G4cout << "A, Z, momentum, theExcitationEnergy"<<
           A<<" "<<Z<<" "<<aMomentum<<" "<<theExcitationEnergy<<G4endl;
      G4Exception( "G4Fragment::G4Fragment Excitation Energy < 0.0!" );
    }
}


// This constructor is for initialize photons
G4Fragment::G4Fragment(const G4LorentzVector aMomentum, G4ParticleDefinition * aParticleDefinition) :
  theA(0),
  theZ(0),
  theMomentum(aMomentum),
  numberOfParticles(0),
  numberOfHoles(0),
  numberOfCharged(0),
  theParticleDefinition(aParticleDefinition),
  theCreationTime(0.0)
#ifdef pctest
  ,theCreatorModel("No name")
#endif
{
  theExcitationEnergy = CalculateExcitationEnergy(aMomentum);
  theAngularMomentum = IsotropicRandom3Vector();
}



const G4Fragment & G4Fragment::operator=(const G4Fragment &right)
{
  if (this != &right) {
    theA = right.theA;
    theZ = right.theZ;
    theExcitationEnergy = right.theExcitationEnergy;
    theMomentum  = right.theMomentum;
    theAngularMomentum = right.theAngularMomentum;
    numberOfParticles = right.numberOfParticles;
    numberOfHoles = right.numberOfHoles;
    numberOfCharged = right.numberOfCharged;
    theParticleDefinition = right.theParticleDefinition;
    theCreationTime = right.theCreationTime;
#ifdef pctest
    theCreatorModel = right.theCreatorModel;
#endif
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


G4std::ostream& operator << (G4std::ostream &out, const G4Fragment *theFragment)
{
#ifdef G4USE_STD_NAMESPACE
  G4std::ios::fmtflags old_floatfield = out.flags();
  out.setf(G4std::ios::floatfield);
#else
  long old_floatfield = out.setf(0,G4std::ios::floatfield);
#endif

  out 
    << "Fragment: A = " << G4std::setprecision(3) << theFragment->theA 
    << ", Z = " << G4std::setprecision(3) << theFragment->theZ ;
  out.setf(G4std::ios::scientific,G4std::ios::floatfield);
  out
    << ", U = " << theFragment->GetExcitationEnergy()/MeV 
    << " MeV" << G4endl
    << "          P = (" 
    << theFragment->theMomentum.x()/MeV << ","
    << theFragment->theMomentum.y()/MeV << ","
    << theFragment->theMomentum.z()/MeV 
    << ") MeV   E = " 
    << theFragment->theMomentum.t()/MeV << " MeV";

  // What about Angular momentum???

  if (theFragment->GetNumberOfExcitons() != 0) {
    out << G4endl;
    out << "          " 
	<< "#Particles = " << theFragment->numberOfParticles 
	<< ", #Holes = "   << theFragment->numberOfHoles
	<< ", #Charged = " << theFragment->numberOfCharged;
  }
  out.setf(old_floatfield,G4std::ios::floatfield);

  return out;
    
}

G4std::ostream& operator << (G4std::ostream &out, const G4Fragment &theFragment)
{
  out << &theFragment;
  return out; 
}



G4double G4Fragment::CalculateExcitationEnergy(const G4LorentzVector value) const
{
	G4double theMaxGroundStateMass = theZ*G4Proton::Proton()->GetPDGMass()+
	                       (theA-theZ)*G4Neutron::Neutron()->GetPDGMass();
	G4double U = value.m() - G4std::min(theMaxGroundStateMass, GetGroundStateMass());
	if( U < 0.0 )
		if( U > -10.0 * eV )
			U = 0.0;
		else 
		{
			G4cerr << "G4Fragment::G4Fragment Excitation Energy ="
			       <<U << " for A = "<<theA<<" and Z= "<<theZ<<G4endl;
			U=0.0;
		}
	return U;
}

G4ThreeVector G4Fragment::IsotropicRandom3Vector(const G4double Magnitude) const
  // Create a unit vector with a random direction isotropically distributed
{

  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*cos(Phi)*SinTheta,
                       Magnitude*sin(Phi)*SinTheta,
                       Magnitude*CosTheta);

  return Vector;
		
}
