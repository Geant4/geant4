// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidate.cc,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCandidate ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Candidates used by CHIPS Model
// ------------------------------------------------------------------

#include "G4QCandidate.hh"

//G4QCandidate::G4QCandidate(G4int PDGcode, G4Nucleus* aNucleus = NULL)
G4QCandidate::G4QCandidate(G4int PDGcode)
{
  PDGencoding = PDGcode;
  theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
  if (theDefinition == NULL)
  {
	cout << "Encoding = " << PDGencoding << G4endl;
	G4Exception("G4QCandidate::GetDefinition(): Encoding not in particle table");
  }
  theMass = theDefinition->GetPDGMass();
  valQ.SetU(theDefinition->GetQuarkContent(2));
  valQ.SetD(theDefinition->GetQuarkContent(1));
  valQ.SetS(theDefinition->GetQuarkContent(3));
  valQ.SetAU(theDefinition->GetAntiQuarkContent(2));
  valQ.SetAD(theDefinition->GetAntiQuarkContent(1));
  valQ.SetAS(theDefinition->GetAntiQuarkContent(3));
  relativeProbability = 0.;
  integralProbability = 0.;
  // Now it's only for Vacuum fragments (Nuclear fragments are initialized differently)
  //if (aNucleus==NULL)
     theMomentum.setE(theMass);
  //else
  //{
	 //G4double z=aNucleus->GetZ();
	 //G4double n=aNucleus->GetN();
	 //cout << "Nucleus M(Z=" << z << ", N=" << n << ")=" << aNucleus->AtomicMass(z+n,z) << G4endl;
	 //G4Exception("G4QCandidate: Initialization for nuclear environment is not implemented");
  //}
}

G4QCandidate::G4QCandidate(const G4QCandidate &right)
{
  theMomentum         = right.theMomentum;
  PDGencoding         = right.PDGencoding;
  theDefinition       = right.theDefinition;
  theMass             = right.theMass;
  valQ                = right.valQ;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
}

const G4QCandidate & G4QCandidate::operator=(const G4QCandidate &right)
{
  theMomentum         = right.theMomentum;
  PDGencoding         = right.PDGencoding;
  theDefinition       = right.theDefinition;
  theMass             = right.theMass;
  valQ                = right.valQ;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
		
  return *this;
}

G4QCandidate::~G4QCandidate() {}

