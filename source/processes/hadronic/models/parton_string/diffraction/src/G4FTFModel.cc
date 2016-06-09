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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4FTFModel.cc,v 1.5 2005/06/04 13:47:01 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
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

G4FTFModel::G4FTFModel(G4DiffractiveExcitation * anExcitation) 
: 
theExcitation(anExcitation)
{
	G4VPartonStringModel::SetThisPointer(this);
}



G4FTFModel::~G4FTFModel()
{}


const G4FTFModel & G4FTFModel::operator=(const G4FTFModel &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFModel::operator= is not meant to be accessed ");
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

struct DeleteVSplitableHadron { void operator()(G4VSplitableHadron * aH){delete aH;} };

G4ExcitedStringVector * G4FTFModel::BuildStrings()
{	

// Loop over all collisions; find all primaries, and all target ( targets may 
//  be duplicate in the List ( to unique G4VSplitableHadrons)

	G4ExcitedStringVector * strings;
	strings = new G4ExcitedStringVector();
	
	std::vector<G4VSplitableHadron *> primaries;
	std::vector<G4VSplitableHadron *> targets;
	
	theParticipants.StartLoop();    // restart a loop 
	while ( theParticipants.Next() ) 
	{
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                 //  do not allow for duplicates ...
	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(), interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());
		
	    if ( targets.end() == std::find(targets.begin(), targets.end(),interaction.GetTarget()) ) 
	    	targets.push_back(interaction.GetTarget());
	}
	    
	
//	G4cout << "BuildStrings prim/targ " << primaries.entries() << " , " <<
//					     targets.entries() << G4endl;


	unsigned int ahadron;
	for ( ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
	    G4bool isProjectile=true;
	    strings->push_back(theExcitation->String(primaries[ahadron], isProjectile));
	}
	for ( ahadron=0; ahadron < targets.size() ; ahadron++)
	{
	    G4bool isProjectile=false;
	    strings->push_back(theExcitation->String(targets[ahadron], isProjectile));
	}

	std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	primaries.clear();
	std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
	targets.clear();
	
	return strings;
}

G4bool G4FTFModel::ExciteParticipants()
{
//	G4cout << "G4FTFModel::ExciteParticipants starting " << G4endl;
	
	while (theParticipants.Next())
	{
//	   G4cout << "next Collision " << G4endl;
	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
	   
//	   G4cout << " soft colls : " << collision.GetNumberOfSoftCollisions() << G4endl;

	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
	   if ( ! theExcitation->ExciteParticipants(projectile, target) ) 
	   {
//           give up, clean up
		std::vector<G4VSplitableHadron *> primaries;
		std::vector<G4VSplitableHadron *> targets;
		theParticipants.StartLoop();    // restart a loop 
		while ( theParticipants.Next() ) 
		{
		    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                	 //  do not allow for duplicates ...
		    if ( primaries.end() == std::find(primaries.begin(), primaries.end(), interaction.GetProjectile()) )
	    		primaries.push_back(interaction.GetProjectile());

		    if ( targets.end() == std::find(targets.begin(), targets.end(),interaction.GetTarget()) ) 
	    		targets.push_back(interaction.GetTarget());
		}
		std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
		primaries.clear();
		std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
		targets.clear();

	   	return false;
	   }

	}
	return true;
}


