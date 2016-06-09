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
// $Id: G4FTFModel.cc,v 1.13 2008/12/09 10:40:52 vuzhinsk Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------

#include "G4FTFModel.hh"
#include "G4FTFParameters.hh"                            // Uzhi 29.03.08
#include "G4FTFParticipants.hh"
#include "G4InteractionContent.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include <utility>                                        // Uzhi 29.03.08

// Class G4FTFModel 

G4FTFModel::G4FTFModel():theExcitation(new G4DiffractiveExcitation()),
                         theElastic(new G4ElasticHNScattering()) // Uzhi 29.03.08
{
	G4VPartonStringModel::SetThisPointer(this);
        theParameters=0;                                         // Uzhi 9.12.08
}

/*
G4FTFModel::G4FTFModel(G4double , G4double , G4double ):theExcitation(new // Uzhi 9.12.08 G4DiffractiveExcitation())
{
	G4VPartonStringModel::SetThisPointer(this);
}

G4FTFModel::G4FTFModel(G4DiffractiveExcitation * anExcitation) 
: 
theExcitation(anExcitation)
{
	G4VPartonStringModel::SetThisPointer(this);
}
*/


G4FTFModel::~G4FTFModel()
{
   if( theParameters != 0 ) delete theParameters;             // Uzhi 5.12.08
// Because FTF model can be called for various particles
// theParameters must be erased at the end of each call.
// Thus the delete is olso in G4FTFModel::GetStrings() method
   if( theExcitation != 0 ) delete theExcitation;             // Uzhi 5.12.08
   if( theElastic    != 0 ) delete theElastic;                // Uzhi 5.12.08
}


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

// ------------------------------------------------------------
void G4FTFModel::Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile)
{
	theProjectile = aProjectile;  
//G4cout<<"G4FTFModel::Init "<<aNucleus.GetN()<<" "<<aNucleus.GetZ()<<G4endl;
	theParticipants.Init(aNucleus.GetN(),aNucleus.GetZ()); 
// Uzhi N-mass number Z-charge ------------------------- Uzhi 29.03.08

// --- cms energy

        G4double s = sqr( theProjectile.GetMass() ) +
                     sqr( G4Proton::Proton()->GetPDGMass() ) +
                     2*theProjectile.GetTotalEnergy()*G4Proton::Proton()->GetPDGMass();
/*
G4cout << " primary Total E (GeV): " << theProjectile.GetTotalEnergy()/GeV << G4endl;
G4cout << " primary Mass    (GeV): " << theProjectile.GetMass() /GeV << G4endl;
G4cout << "cms std::sqrt(s) (GeV) = " << std::sqrt(s) / GeV << G4endl;
*/

      if( theParameters != 0 ) delete theParameters;                    // Uzhi 9.12.08
      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
                                          aNucleus.GetN(),aNucleus.GetZ(),
                                          s);// ------------------------- Uzhi 19.04.08
//theParameters->SetProbabilityOfElasticScatt(0.); // To turn on/off (1/0) elastic scattering
}

// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::GetStrings()
{
//G4cout<<"theParticipants.GetList"<<G4endl;
	theParticipants.GetList(theProjectile,theParameters);
//G4cout<<"ExciteParticipants()"<<G4endl;
	if (! ExciteParticipants()) return NULL;;
//G4cout<<"theStrings = BuildStrings()"<<G4endl;
	G4ExcitedStringVector * theStrings = BuildStrings();
//G4cout<<"Return to theStrings "<<G4endl;
        if( theParameters != 0 )                              // Uzhi 9.12.08 
        {                                                     // Uzhi 9.12.08
          delete theParameters;                               // Uzhi 9.12.08
          theParameters=0;                                    // Uzhi 9.12.08
        }                                                     // Uzhi 9.12.08
	return theStrings;
}

// ------------------------------------------------------------
struct DeleteVSplitableHadron { void operator()(G4VSplitableHadron * aH){delete aH;} };

// ------------------------------------------------------------
G4bool G4FTFModel::ExciteParticipants()
{
/*    // Uzhi 29.03.08                     For elastic Scatt.
G4cout<<"  In ExciteParticipants() "<<theParticipants.theInteractions.size()<<G4endl;
G4cout<<" test Params Tot "<<theParameters->GetTotalCrossSection()<<G4endl;
G4cout<<" test Params Ela "<<theParameters->GetElasticCrossSection()<<G4endl;
	
G4int counter=0;
*/   // Uzhi 29.03.08



	while (theParticipants.Next())
	{	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
/*
counter++;
G4cout<<" Inter # "<<counter<<G4endl;
*/
	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();

//   // Uzhi 29.03.08
           G4bool Successfull;
           if(G4UniformRand()< theParameters->GetProbabilityOfElasticScatt())
           {
//G4cout<<"Elastic"<<G4endl;
            Successfull=theElastic->ElasticScattering(projectile, target, theParameters);
           }
           else
           {
//G4cout<<"Inelastic"<<G4endl;
            Successfull=theExcitation->ExciteParticipants(projectile, target, theParameters);
           }
//           if(!Successfull)
// // Uzhi 29.03.08

//	   if ( ! theExcitation->ExciteParticipants(projectile, target) ) 
           if(!Successfull)
	   {
//           give up, clean up
		std::vector<G4VSplitableHadron *> primaries;
		std::vector<G4VSplitableHadron *> targets;
		theParticipants.StartLoop();    // restart a loop 
		while ( theParticipants.Next() ) 
		{
		    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                	 //  do not allow for duplicates ...
		    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                      interaction.GetProjectile()) )
	    		primaries.push_back(interaction.GetProjectile());

		    if ( targets.end()   == std::find(targets.begin(), targets.end(),
                                                      interaction.GetTarget()) ) 
	    		targets.push_back(interaction.GetTarget());
		}
		std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
		primaries.clear();
		std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
		targets.clear();

	   	return false;
	   }  // End of the loop Uzhi
        }
	return true;
}
// ------------------------------------------------------------
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
	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());
		
	    if ( targets.end()   == std::find(targets.begin(), targets.end(),
                                                interaction.GetTarget()) ) 
	    	targets.push_back(interaction.GetTarget());
	}
	    
	
//	G4cout << "BuildStrings prim/targ " << primaries.size() << " , " <<
//					       targets.size() << G4endl;

	unsigned int ahadron;
// Only for hA-interactions Uzhi -------------------------------------
	for ( ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
//G4ThreeVector aPosition=primaries[ahadron]->GetPosition();
//G4cout<<"Proj Build "<<aPosition<<" "<<primaries[ahadron]->GetTimeOfCreation()<<G4endl;
	    G4bool isProjectile=true;
	    strings->push_back(theExcitation->String(primaries[ahadron], isProjectile));
	}

	for ( ahadron=0; ahadron < targets.size() ; ahadron++)
	{
//G4ThreeVector aPosition=targets[ahadron]->GetPosition();
//G4cout<<"Targ Build "<<aPosition<<" "<<targets[ahadron]->GetTimeOfCreation()<<G4endl;
	    G4bool isProjectile=false;
	    strings->push_back(theExcitation->String(targets[ahadron], isProjectile));
	}

	std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	primaries.clear();
	std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
	targets.clear();
	
	return strings;
}
// ------------------------------------------------------------
