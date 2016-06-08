// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPartonStringModel.cc,v 1.7 2000/08/02 08:13:35 hpw Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
//// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VPartonStringModel ----------------
//             by Gunter Folger, May 1998.
//      abstract class for all Parton String Models
// ------------------------------------------------------------


#include "G4VPartonStringModel.hh"
#include "G4ios.hh"
#include "G4ShortLivedConstructor.hh"


G4VPartonStringModel::G4VPartonStringModel()
{
//  Make shure Shotrylived partyicles are constructed.
	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
}

G4VPartonStringModel::G4VPartonStringModel(const G4VPartonStringModel &right)
{
}


G4VPartonStringModel::~G4VPartonStringModel()
{
}


const G4VPartonStringModel & G4VPartonStringModel::operator=(const G4VPartonStringModel &right)
{
  G4Exception("G4VPartonStringModel::operator= meant to not be accessable");
  return *this;
}


int G4VPartonStringModel::operator==(const G4VPartonStringModel &right) const
{
 return 0;
}

int G4VPartonStringModel::operator!=(const G4VPartonStringModel &right) const
{
  return 1;
}

G4KineticTrackVector * G4VPartonStringModel::Scatter(const G4Nucleus &theNucleus, 
                                                const G4DynamicParticle &aPrimary)
{  
  G4ExcitedStringVector * strings = NULL;

  G4DynamicParticle thePrimary=aPrimary;
  
  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=thePrimary.Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  thePrimary.Set4Momentum(toZ*Ptmp);
  G4LorentzRotation toLab(toZ.inverse());

  G4int attempts = 0, maxAttempts=20;
  while ( strings  == NULL )
  {
  	if (attempts++ > maxAttempts ) 
  	{
		G4Exception("G4VPartonStringModel::Scatter(): fails to generate strings");
  	}

	theThis->Init(theNucleus,thePrimary);
  	strings = GetStrings();
  }
  
  G4KineticTrackVector * theResult = 0;
  for ( G4int astring=0; astring < strings->entries(); astring++)
  {
//    rotate string to lab frame, models have it aligned to z
  	strings->at(astring)->LorentzRotate(toLab);
  }
  
  theResult = stringFragmentationModel->FragmentStrings(strings);
  strings->clearAndDestroy();
  delete strings;

  return theResult;
}

