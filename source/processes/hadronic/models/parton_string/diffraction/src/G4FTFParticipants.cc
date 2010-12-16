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
// $Id: G4FTFParticipants.cc,v 1.17 2010-09-20 15:50:46 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
//  Changed in a part by V. Uzhinsky in oder to put in correcpondence
//        with original FRITIOF mode. November - December 2006.
// ------------------------------------------------------------

#include "G4FTFParameters.hh"                            // Uzhi 29.03.08
#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4VSplitableHadron.hh"
#include "Randomize.hh"
#include <utility>                                        // Uzhi 29.03.08

// Class G4FTFParticipants 

G4FTFParticipants::G4FTFParticipants() 
{
}

G4FTFParticipants::G4FTFParticipants(const G4FTFParticipants &): G4VParticipants()
{
}


G4FTFParticipants::~G4FTFParticipants()
{
}


//const G4FTFParticipants & G4FTFParticipants::operator=(const G4FTFParticipants &right)
//{}


//int G4FTFParticipants::operator==(const G4FTFParticipants &right) const
//{}

//int G4FTFParticipants::operator!=(const G4FTFParticipants &right) const
//{}

void G4FTFParticipants::GetList(const G4ReactionProduct  &thePrimary,
                                      G4FTFParameters    *theParameters) 
{
    
    StartLoop();  // reset Loop over Interactions

    for(unsigned int i=0; i<theInteractions.size(); i++) delete theInteractions[i];
    theInteractions.clear();

    G4double deltaxy=2 * fermi;                       // Extra nuclear radius

    G4VSplitableHadron * primarySplitable=new G4DiffractiveSplitableHadron(thePrimary);

    G4double xyradius;                          
    xyradius =theNucleus->GetOuterRadius() + deltaxy; // Impact parameter sampling
                                                      // radius
//    G4bool nucleusNeedsShift = true;                // Uzhi 20 July 2009
    
    while ( theInteractions.size() == 0 )
    {
	std::pair<G4double, G4double> theImpactParameter;
	theImpactParameter = theNucleus->ChooseImpactXandY(xyradius);
	G4double impactX = theImpactParameter.first; 
	G4double impactY = theImpactParameter.second;

        G4ThreeVector thePosition(impactX, impactY, -DBL_MAX);     
        primarySplitable->SetPosition(thePosition);                

	theNucleus->StartLoop();
	G4Nucleon * nucleon;

	while ( (nucleon=theNucleus->GetNextNucleon()) ) 
	{
    	   G4double impact2= sqr(impactX - nucleon->GetPosition().x()) +
                             sqr(impactY - nucleon->GetPosition().y());

	   if ( theParameters->GetProbabilityOfInteraction(impact2/fermi/fermi) 
		> G4UniformRand() )
	   {
                primarySplitable->SetStatus(1);        // It takes part in the interaction

		G4VSplitableHadron * targetSplitable=0;
	   	if ( ! nucleon->AreYouHit() )
	   	{
	   	    targetSplitable= new G4DiffractiveSplitableHadron(*nucleon);
	   	    nucleon->Hit(targetSplitable);
	   	    nucleon->SetBindingEnergy(3.*nucleon->GetBindingEnergy()); 
                    targetSplitable->SetStatus(1);     // It takes part in the interaction
	   	}
	   	G4InteractionContent * aInteraction = 
                                       new G4InteractionContent(primarySplitable);
		aInteraction->SetTarget(targetSplitable);
                aInteraction->SetTargetNucleon(nucleon);     // Uzhi 16.07.09
		theInteractions.push_back(aInteraction);
	   }
	}    

//	G4cout << "Number of Hit nucleons " << theInteractions.size()<<G4endl; //  entries() 
//		<< "\t" << impactX/fermi << "\t"<<impactY/fermi
//		<< "\t" << std::sqrt(sqr(impactX)+sqr(impactY))/fermi <<G4endl;

    }
   
}


// Implementation (private) methods
