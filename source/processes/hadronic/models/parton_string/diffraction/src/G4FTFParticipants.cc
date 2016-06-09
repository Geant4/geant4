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
// $Id$
// GEANT4 tag $Name:  $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
//  Changed in a part by V. Uzhinsky in oder to put in correcpondence
//        with original FRITIOF mode. November - December 2006.
//  Ajusted for (anti) nucleus - nucleus interactions by V. Uzhinsky.
//                    (February 2011)
// ------------------------------------------------------------

#include <utility>
#include <vector>
#include <algorithm>

#include "G4FTFParticipants.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4FTFParameters.hh"                            // Uzhi 29.03.08
#include "G4DiffractiveSplitableHadron.hh"
#include "G4VSplitableHadron.hh"

G4FTFParticipants::G4FTFParticipants() :
  theProjectileNucleus(0),
  currentInteraction(-1)
{
}

G4FTFParticipants::G4FTFParticipants(const G4FTFParticipants &): G4VParticipants()
  , theProjectileNucleus(0), currentInteraction(-1)   //A.R. 14-Aug-2012 Coverity fix.
{
	G4Exception("G4FTFParticipants::G4FTFParticipants()","HAD_FTF_001",
	        FatalException," Must not use copy ctor()");
}


G4FTFParticipants::~G4FTFParticipants()
{
	if ( theProjectileNucleus != NULL ) delete theProjectileNucleus;
}

//-------------------------------------------------------------------------

void G4FTFParticipants::SetProjectileNucleus(G4V3DNucleus * aNucleus)
{
  if (theProjectileNucleus) delete theProjectileNucleus;
  
  theProjectileNucleus = aNucleus;
}

G4V3DNucleus * G4FTFParticipants::GetProjectileNucleus()
{
  return theProjectileNucleus;
}

void G4FTFParticipants::InitProjectileNucleus(G4int theA, G4int theZ)
{
	if ( theProjectileNucleus == NULL ) theProjectileNucleus = new G4Fancy3DNucleus();
	theProjectileNucleus->Init(theA, theZ);
        theProjectileNucleus->SortNucleonsDecZ();
}
//-------------------------------------------------------------------------
void G4FTFParticipants::GetList(const G4ReactionProduct  &thePrimary,
                                      G4FTFParameters    *theParameters) 
{ 
//G4cout<<"Participants::GetList"<<G4endl;
//G4cout<<"thePrimary "<<thePrimary.GetMomentum()<<G4endl;
    StartLoop();  // reset Loop over Interactions

    for(unsigned int i=0; i<theInteractions.size(); i++) delete theInteractions[i];
    theInteractions.clear();

    G4double deltaxy=2 * fermi;                       // Extra nuclear radius
//G4cout<<"theProjectileNucleus "<<theProjectileNucleus<<G4endl;
    if(theProjectileNucleus == 0)
    { // Hadron-nucleus or anti-baryon-nucleus interactions
//G4cout<<"Hadron-nucleus or anti-baryon-nucleus interactions"<<G4endl;

     G4double impactX(0.), impactY(0.);

     G4VSplitableHadron * primarySplitable=new G4DiffractiveSplitableHadron(thePrimary);
//G4cout<<"Prim in Part "<<primarySplitable->Get4Momentum()<<G4endl;
     G4double xyradius;                          
     xyradius =theNucleus->GetOuterRadius() + deltaxy; // Impact parameter sampling
                                                    
//    G4bool nucleusNeedsShift = true;                // Uzhi 20 July 2009
    
     do 
     {
	 std::pair<G4double, G4double> theImpactParameter;
	 theImpactParameter = theNucleus->ChooseImpactXandY(xyradius);
	 impactX = theImpactParameter.first; 
	 impactY = theImpactParameter.second;

         G4ThreeVector thePosition(impactX, impactY, -DBL_MAX);     
         primarySplitable->SetPosition(thePosition);                

	 theNucleus->StartLoop();
	 G4Nucleon * nucleon;

G4int TrN(0);
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
//G4cout<<" Part nucl "<<TrN<<" "<<nucleon->Get4Momentum()<<G4endl;
//G4cout<<" Part nucl "<<G4endl;
                    targetSplitable->SetStatus(1);     // It takes part in the interaction
	   	}
	   	G4InteractionContent * aInteraction = 
                                       new G4InteractionContent(primarySplitable);
		aInteraction->SetTarget(targetSplitable);
                aInteraction->SetTargetNucleon(nucleon);     // Uzhi 16.07.09
                aInteraction->SetStatus(1);                  // Uzhi Feb26
		theInteractions.push_back(aInteraction);
	   }
TrN++;
	 } 
     } while ( theInteractions.size() == 0 );

     //if ( theInteractions.size() == 0 ) delete primarySplitable; //A.R. 14-Aug-2012 Coverity fix

//	G4cout << "Number of Hit nucleons " << theInteractions.size()
//		<< "\t" << impactX/fermi << "\t"<<impactY/fermi
//		<< "\t" << std::sqrt(sqr(impactX)+sqr(impactY))/fermi <<G4endl;

     return;
    }       // end of if(theProjectileNucleus == 0)

//-------------------------------------------------------------------
//                Projectile and target are nuclei
//-------------------------------------------------------------------
//VU    G4VSplitableHadron * primarySplitable=new G4DiffractiveSplitableHadron(thePrimary);
//G4cout<<"Prim in Part "<<primarySplitable->Get4Momentum()<<G4endl;
//G4cout<<"Projectile and target are nuclei"<<G4endl;
//G4cout<<thePrimary.GetMomentum()<<G4endl;
//G4cout<<"Part Pr Tr "<<theProjectileNucleus<<" "<<theNucleus<<G4endl;


    G4double xyradius;                          
    xyradius =theProjectileNucleus->GetOuterRadius() +  // Impact parameter sampling
                        theNucleus->GetOuterRadius() + deltaxy;

    G4double impactX(0.), impactY(0.);

    do
    {
//G4cout<<"New interaction list"<<G4endl;
	 std::pair<G4double, G4double> theImpactParameter;
	 theImpactParameter = theNucleus->ChooseImpactXandY(xyradius);
	 impactX = theImpactParameter.first; 
	 impactY = theImpactParameter.second;
//G4cout<<"B "<<std::sqrt(sqr(impactX)+sqr(impactY))/fermi<<G4endl;

         G4ThreeVector thePosition(impactX, impactY, -DBL_MAX);     
//VU         primarySplitable->SetPosition(thePosition);                

	 theProjectileNucleus->StartLoop();
	 G4Nucleon * ProjectileNucleon;
G4int PrNuclN(0);

	 while ( (ProjectileNucleon=theProjectileNucleus->GetNextNucleon()) ) 
	 {
           G4VSplitableHadron * ProjectileSplitable=0;
//G4cout<<G4endl<<"Prj N mom "<<ProjectileNucleon->Get4Momentum()<<"-------------"<<G4endl;
           theNucleus->StartLoop();
           G4Nucleon * TargetNucleon;

G4int TrNuclN(0);
           while ( (TargetNucleon=theNucleus->GetNextNucleon()) )
           {
//G4cout<<"Trg N mom "<<TargetNucleon->Get4Momentum()<<G4endl;
    	    G4double impact2=
            sqr(impactX+ProjectileNucleon->GetPosition().x()-TargetNucleon->GetPosition().x())+
            sqr(impactY+ProjectileNucleon->GetPosition().y()-TargetNucleon->GetPosition().y());

            G4VSplitableHadron * TargetSplitable=0;

	    if ( theParameters->GetProbabilityOfInteraction(impact2/fermi/fermi) 
		 > G4UniformRand() )
	    { // An Interaction has happend!
//G4cout<<"An Interaction has happend"<<G4endl;
//G4cout<<"PrN TrN "<<PrNuclN<<" "<<TrNuclN<<" "<<ProjectileNucleon->GetPosition().z()/fermi<<" "<<TargetNucleon->GetPosition().z()/fermi<<" "<<ProjectileNucleon->GetPosition().z()/fermi + TargetNucleon->GetPosition().z()/fermi <<G4endl;

             if ( ! ProjectileNucleon->AreYouHit() )
             { // Projectile nucleon was not involved until now.
              ProjectileSplitable= new G4DiffractiveSplitableHadron(*ProjectileNucleon);
              ProjectileNucleon->Hit(ProjectileSplitable);
              ProjectileNucleon->SetBindingEnergy(3.*ProjectileNucleon->GetBindingEnergy());
              ProjectileSplitable->SetStatus(1);     // It takes part in the interaction
             }
             else
             {  // Projectile nucleon was involved before.
              ProjectileSplitable=ProjectileNucleon->GetSplitableHadron();
             } // End of if ( ! Projectileucleon->AreYouHit() )

             if ( ! TargetNucleon->AreYouHit() )
             {  // Target nucleon was not involved until now
              TargetSplitable= new G4DiffractiveSplitableHadron(*TargetNucleon);
              TargetNucleon->Hit(TargetSplitable);
              TargetNucleon->SetBindingEnergy(3.*ProjectileNucleon->GetBindingEnergy());
              TargetSplitable->SetStatus(1);     // It takes part in the interaction
             }
             else
             {  // Target nucleon was involved before.
              TargetSplitable=TargetNucleon->GetSplitableHadron();
             } // End of if ( ! TargetNeucleon->AreYouHit() )

             G4InteractionContent * anInteraction = 
                                   new G4InteractionContent(ProjectileSplitable);
             anInteraction->SetTarget(TargetSplitable);
             anInteraction->SetTargetNucleon(TargetNucleon);
             anInteraction->SetStatus(1);                      // Uzhi Feb26
//             anInteraction->SetInteractionTime(ProjectileNucleon->GetPosition().z()+
//                                                   TargetNucleon->GetPosition().z());
//G4cout<<"Z's pr tr "<<ProjectileNucleon->GetPosition().z()/fermi<<" "<<TargetNucleon->GetPosition().z()/fermi<<" "<<ProjectileNucleon->GetPosition().z()/fermi + TargetNucleon->GetPosition().z()/fermi <<G4endl;
             theInteractions.push_back(anInteraction);
//G4cout<<"Ppr tr "<<ProjectileSplitable<<" "<<TargetSplitable<<G4endl;
            } // End of An Interaction has happend!
TrNuclN++;
           } // End of while ( (TargetNucleon=theNucleus->GetNextNucleon()) )
PrNuclN++;
	 } // End of   while ( (ProjectileNucleon=theProjectileNucleus->GetNextNucleon()) )
    }  while ( theInteractions.size() == 0 );  // end of while ( theInteractions.size() == 0 )

//std::sort(theInteractions.begin(),theInteractions.end()); // ????

//	G4cout << "Number of primary collisions " << theInteractions.size() 
//		<< "\t" << impactX/fermi << "\t"<<impactY/fermi
//		<< "\t" << std::sqrt(sqr(impactX)+sqr(impactY))/fermi <<G4endl;
//G4int Uzhi; G4cin >> Uzhi;
    return;
}
//--------------------------------------------------------------

// Implementation (private) methods
