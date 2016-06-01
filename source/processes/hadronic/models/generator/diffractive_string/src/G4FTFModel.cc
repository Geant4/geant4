// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FTFModel.cc,v 1.5 1998/12/09 07:41:26 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------

#include "G4FTFModel.hh"
#include "G4FTFParticipants.hh"
#include "G4InteractionContent.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

// Class G4FTFModel 

G4FTFModel::G4FTFModel(G4double sigmaPt, G4double minExtraMass,G4double x0Mass) 
: 
theExcitation(new G4DiffractiveExcitation(sigmaPt,minExtraMass,x0Mass))
{
	G4VPartonStringModel::SetThisPointer(this);
}



G4FTFModel::~G4FTFModel()
{}


const G4FTFModel & G4FTFModel::operator=(const G4FTFModel &right)
{
	G4Exception("G4FTFModel::operator= is not meant to be accessed ");
	return *this;
}


int G4FTFModel::operator==(const G4FTFModel &right) const
{
	return this==&right;
}

int G4FTFModel::operator!=(const G4FTFModel &right) const
{
	return this!=&right;
}

void G4FTFModel::Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile)
{
	theParticipants.Init(aNucleus.GetN(),aNucleus.GetZ());
	theProjectile = aProjectile;  
}

G4ExcitedStringVector * G4FTFModel::GetStrings()
{

	theParticipants.BuildInteractions(theProjectile);
	
	if (! ExciteParticipants()) return NULL;;

	G4ExcitedStringVector * theStrings = BuildStrings();

	return theStrings;
}

G4ExcitedStringVector * G4FTFModel::BuildStrings()
{	

// Loop over all collisions; find all primaries, and all target ( targets may 
//  be duplicate in the List ( to unique G4VSplitableHadrons)

	G4ExcitedStringVector * strings;
	strings = new G4ExcitedStringVector();
	
	RWTPtrOrderedVector<G4VSplitableHadron> primaries;
	RWTPtrOrderedVector<G4VSplitableHadron> targets;
	
	theParticipants.StartLoop();    // restart a loop 
	while ( theParticipants.Next() ) 
	{
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                 //  do not allow for duplicates ...
	    if ( ! primaries.contains(interaction.GetProjectile()) )
	    	primaries.insert(interaction.GetProjectile());
	    if ( ! targets.contains(interaction.GetTarget()) ) 
	    	targets.insert(interaction.GetTarget());
	}
	    
	
//	G4cout << "BuildStrings prim/targ " << primaries.entries() << " , " <<
//					     targets.entries() << endl;


	G4int ahadron;
	for ( ahadron=0; ahadron < primaries.entries() ; ahadron++)
	{
	    G4bool isProjectile=true;
	    strings->insert(theExcitation->String(primaries[ahadron], isProjectile));
	}
	for ( ahadron=0; ahadron < targets.entries() ; ahadron++)
	{
	    G4bool isProjectile=false;
	    strings->insert(theExcitation->String(targets[ahadron], isProjectile));
	}

	primaries.clearAndDestroy();
	targets.clearAndDestroy();
	
	return strings;
}

G4bool G4FTFModel::ExciteParticipants()
{
//	G4cout << "G4FTFModel::ExciteParticipants starting " << endl;
	
	while (theParticipants.Next())
	{
//	   G4cout << "next Collision " << endl;
	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
	   
//	   G4cout << " soft colls : " << collision.GetNumberOfSoftCollisions() << endl;

	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
	   if ( ! theExcitation->ExciteParticipants(projectile, target) ) 
	   {
	   	return false;
	   }

	}
	return true;
}


