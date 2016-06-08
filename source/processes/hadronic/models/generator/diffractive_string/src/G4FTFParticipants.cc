// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FTFParticipants.cc,v 1.3 2000/08/19 16:23:48 hpw Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
// ------------------------------------------------------------

#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4VSplitableHadron.hh"
#include "G4PomeronCrossSection.hh"
#include "Randomize.hh"
#include "G4Pair.hh"


// Class G4FTFParticipants 



G4FTFParticipants::G4FTFParticipants() 
{
}

G4FTFParticipants::G4FTFParticipants(const G4FTFParticipants &right)
{
}


G4FTFParticipants::~G4FTFParticipants()
{
// G4cout << "G4FTFParticipants::~G4FTFParticipants() called" << G4endl;
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

    theInteractions.clearAndDestroy();

// --- cms energy

    G4double s = sqr( thePrimary.GetMass() ) +
		 sqr( G4Proton::Proton()->GetPDGMass() ) +
		 2*thePrimary.GetTotalEnergy()*G4Proton::Proton()->GetPDGMass();

//    G4cout << " primary Total E (GeV): " << thePrimary.GetTotalEnergy()/GeV << G4endl;
//    G4cout << " primary Mass    (GeV): " << thePrimary.GetMass() /GeV << G4endl;
//    G4cout << "cms sqrt(s) (GeV) = " << sqrt(s) / GeV << G4endl;

    G4PomeronCrossSection theCrossSection(thePrimary.GetDefinition());
    
    

    G4double deltaxy=2 * fermi; 

    G4VSplitableHadron * primarySplitable=new G4DiffractiveSplitableHadron(thePrimary);

    G4double xyradius;
    xyradius =theNucleus->GetOuterRadius() + deltaxy;

//    G4cout <<"  G4FTFParticipants::StartLoop: xyradius " << xyradius << G4endl;

    G4bool nucleusNeedsShift = true;
    
    while ( theInteractions.entries() == 0 )
    {
	G4Pair<G4double, G4double> theImpactParameter;
	theImpactParameter = theNucleus->ChooseImpactXandY(xyradius);
	G4double impactX = theImpactParameter.first; 
	G4double impactY = theImpactParameter.second;

//	G4cout << " impctX, impctY " << impactX/fermi << "    "<<impactY/fermi << " fm" << G4endl;

	theNucleus->StartLoop();
	G4Nucleon * nucleon;
	while ( nucleon=theNucleus->GetNextNucleon() )
	{
    	   G4double impact2= sqr(impactX - nucleon->GetPosition().x()) +
    		    sqr(impactY - nucleon->GetPosition().y());
	   if ( theCrossSection.GetInelasticProbability(s,impact2)
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
		theInteractions.insert(aInteraction);
	   }
	}    

//	G4cout << "Number of Hit nucleons " << theInteractions.entries() 
//		<< "\t" << impactX/fermi << "\t"<<impactY/fermi
//		<< "\t" << sqrt(sqr(impactX)+sqr(impactY))/fermi <<G4endl;
	 
    }
   
}


// Implementation (private) methods





