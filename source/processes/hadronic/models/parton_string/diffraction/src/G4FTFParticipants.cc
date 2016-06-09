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
// $Id: G4FTFParticipants.cc,v 1.7 2007/04/24 10:33:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $
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

#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4VSplitableHadron.hh"
//#include "G4PomeronCrossSection.hh"                      // Uzhi
#include "G4FTFCrossSection.hh"                            // Uzhi
#include "Randomize.hh"
#include <utility>


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

void G4FTFParticipants::BuildInteractions(const G4ReactionProduct  &thePrimary)
{
    
    StartLoop();  // reset Loop over Interactions

    for(unsigned int i=0; i<theInteractions.size(); i++) delete theInteractions[i];
    theInteractions.clear();

// --- cms energy

    G4double s = sqr( thePrimary.GetMass() ) +
		 sqr( G4Proton::Proton()->GetPDGMass() ) +
		 2*thePrimary.GetTotalEnergy()*G4Proton::Proton()->GetPDGMass();

//    G4cout << " primary Total E (GeV): " << thePrimary.GetTotalEnergy()/GeV << G4endl;
//    G4cout << " primary Mass    (GeV): " << thePrimary.GetMass() /GeV << G4endl;
//    G4cout << "cms std::sqrt(s) (GeV) = " << std::sqrt(s) / GeV << G4endl;

//    G4PomeronCrossSection theCrossSection(thePrimary.GetDefinition());      // Uzhi
      G4FTFCrossSection theCrossSection(thePrimary.GetDefinition(),s);        // Uzhi
    

    G4double deltaxy=2 * fermi; 

    G4VSplitableHadron * primarySplitable=new G4DiffractiveSplitableHadron(thePrimary);

    G4double xyradius;
    xyradius =theNucleus->GetOuterRadius() + deltaxy;

    G4bool nucleusNeedsShift = true;
    
    while ( theInteractions.size() == 0 )
    {
	std::pair<G4double, G4double> theImpactParameter;
	theImpactParameter = theNucleus->ChooseImpactXandY(xyradius);
	G4double impactX = theImpactParameter.first; 
	G4double impactY = theImpactParameter.second;

	theNucleus->StartLoop();
	G4Nucleon * nucleon;
	while ( (nucleon=theNucleus->GetNextNucleon()) )
	{
    	   G4double impact2= sqr(impactX - nucleon->GetPosition().x()) +
    		    sqr(impactY - nucleon->GetPosition().y());
//	   if ( theCrossSection.GetInelasticProbability(s,impact2)                       // Uzhi
	   if ( theCrossSection.GetInelasticProbability(  impact2/fermi/fermi)           // Uzhi 
		> G4UniformRand() )
	   {
	   	if ( nucleusNeedsShift ) 
	   	{			// on the first hit, shift nucleus 
	   	     nucleusNeedsShift = false;
	   	     theNucleus->DoTranslation(G4ThreeVector(-1*impactX,-1*impactY,0.));
	   	     impactX=0;
	   	     impactY=0;
	   	}
		G4VSplitableHadron * targetSplitable;
	   	if ( (targetSplitable=nucleon->GetSplitableHadron()) == NULL )
	   	{
	   	    targetSplitable= new G4DiffractiveSplitableHadron(*nucleon);
	   	    nucleon->Hit(targetSplitable);
	   	}
	   	G4InteractionContent * aInteraction = new G4InteractionContent(primarySplitable);
		aInteraction->SetTarget(targetSplitable);
		theInteractions.push_back(aInteraction);
	   }
	}    

//	G4cout << "Number of Hit nucleons " << theInteractions.size() //  entries() 
//		<< "\t" << impactX/fermi << "\t"<<impactY/fermi
//		<< "\t" << std::sqrt(sqr(impactX)+sqr(impactY))/fermi <<G4endl;

    }
   
}


// Implementation (private) methods





