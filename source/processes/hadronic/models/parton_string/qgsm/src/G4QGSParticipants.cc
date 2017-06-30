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
#include <utility>

#include "G4QGSParticipants.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4LorentzVector.hh"
#include "G4V3DNucleus.hh" 
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4PhysicalConstants.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//#define debugQGSParticipants
//#define debugPutOnMassShell

// Class G4QGSParticipants 

// Promoting model parameters from local variables class properties
G4ThreadLocal G4int G4QGSParticipants_NPart = 0;

G4QGSParticipants::G4QGSParticipants() : theDiffExcitaton(), //0.7*GeV, 250*MeV, 250*MeV),
		ModelMode(SOFT),
		//nCutMax(7),ThresholdParameter(0.45*GeV),
		nCutMax(7),ThresholdParameter(0.000*GeV),
		QGSMThreshold(3*GeV),theNucleonRadius(1.5*fermi),alpha(-0.5),beta(2.5)
// Are ThresholdParameter and QGSMThreshold needed here?
{
// Parameters setting
  SetCofNuclearDestruction(1.);
  SetR2ofNuclearDestruction( 1.5*fermi*fermi );
  SetDofNuclearDestruction( 0.3 );
  SetPt2ofNuclearDestruction( 0.075*GeV*GeV );
  SetMaxPt2ofNuclearDestruction( 1.0*GeV*GeV );
  SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );

  sigmaPt=0.25*sqr(GeV);   // Uzhi 2 June 2016   
}

G4QGSParticipants::G4QGSParticipants(const G4QGSParticipants &right)
: G4VParticipants(),ModelMode(right.ModelMode), nCutMax(right.nCutMax),
  ThresholdParameter(right.ThresholdParameter), QGSMThreshold(right.QGSMThreshold),
  theNucleonRadius(right.theNucleonRadius)
{
// From FTF Model  Make right Copy
}

G4QGSParticipants::~G4QGSParticipants() {}

void G4QGSParticipants::BuildInteractions(const G4ReactionProduct  &thePrimary) 
{
  theProjectile = thePrimary;

  Regge = new G4Reggeons(theProjectile.GetDefinition()); // Uzhi Oct. 2016

  SetProjectileNucleus( 0 ); // Uzhi theParticipants.SetProjectileNucleus( 0 );

  NumberOfInvolvedNucleonsOfProjectile= 0;
  G4LorentzVector tmp( 0.0, 0.0, 0.0, 0.0 );
  ProjectileResidualMassNumber       = 0;
  ProjectileResidualCharge           = 0;
  ProjectileResidualExcitationEnergy = 0.0;
  ProjectileResidual4Momentum        = tmp;

  NumberOfInvolvedNucleonsOfTarget= 0;
  TargetResidualMassNumber       = theNucleus->GetMassNumber();
  TargetResidualCharge           = theNucleus->GetCharge();
  TargetResidualExcitationEnergy = 0.0;

  theNucleus->StartLoop();
  G4Nucleon* NuclearNucleon;
  while ( ( NuclearNucleon = theNucleus->GetNextNucleon() ) ) {  
    tmp+=NuclearNucleon->Get4Momentum();
  } 
  TargetResidual4Momentum        = tmp;

  if ( std::abs( theProjectile.GetDefinition()->GetBaryonNumber() ) <= 1 ) { 
    // Projectile is a hadron : meson or baryon
    ProjectileResidualMassNumber = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
    ProjectileResidualCharge = G4int( theProjectile.GetDefinition()->GetPDGCharge() );
    ProjectileResidualExcitationEnergy = 0.0;
    ProjectileResidual4Momentum.setVect( theProjectile.GetMomentum() );
    ProjectileResidual4Momentum.setE( theProjectile.GetTotalEnergy() );
  } 

  #ifdef debugQGSParticipants
    G4cout <<G4endl<< "G4QGSParticipants::BuildInteractions---------" << G4endl
           << "thePrimary " << thePrimary.GetDefinition()->GetParticleName()<<" "
           <<ProjectileResidual4Momentum<<G4endl;
    G4cout << "Target A and Z  " << theNucleus->GetMassNumber() <<" "<<theNucleus->GetCharge()<<" "
           << TargetResidual4Momentum<<G4endl;
  #endif

  G4bool Success( true );

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do  // while( (!Success) && ++loopCounter < maxNumberOfLoops );
  {
   const G4int maxNumberOfInternalLoops = 1000;
   G4int internalLoopCounter = 0;
   do  // while( (!Success) && ++internalLoopCounter < maxNumberOfInternalLoops );
   {
    if(std::abs(theProjectile.GetDefinition()->GetPDGEncoding()) < 100.0) 
      {SelectInteractions(theProjectile);}                            // for lepton projectile
    else
      {GetList( theProjectile );  // Get list of participating nucleons  for hadron projectile
    }

    StoreInvolvedNucleon();       // Store participating nucleon

    ReggeonCascade();             // Make reggeon cascading. Involve nucleons in the cascading.

    Success = PutOnMassShell();   // Ascribe momenta to the involved and participating nucleons

    if(!Success) PrepareInitialState( thePrimary );

   } while( (!Success) && ++internalLoopCounter < maxNumberOfInternalLoops );

   if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
     Success = false;
   }
   
   if ( Success ) {
     #ifdef debugQGSParticipants
       G4cout<<"PerformDiffractiveCollisions(); if they happend." <<G4endl;
     #endif

     PerformDiffractiveCollisions();

     #ifdef debugQGSParticipants
       G4cout<<"SplitHadrons();" <<G4endl;
     #endif

     SplitHadrons(); 

     if( theProjectileSplitable && theProjectileSplitable->GetStatus() == 0) {
       #ifdef debugQGSParticipants
         G4cout<<"Perform non-Diffractive Collisions if they happend. Determine Parton Momenta and so on." <<G4endl;
       #endif 
       Success = DeterminePartonMomenta(); 
     }

     if(!Success) PrepareInitialState( thePrimary );
   }
  } while( (!Success) && ++loopCounter < maxNumberOfLoops );

  if ( loopCounter >= maxNumberOfLoops ) {
    Success = false;
    #ifdef debugQGSParticipants
      G4cout<<"NOT Successful ======" <<G4endl;
    #endif
  }

  if ( Success ) {
    CreateStrings();  // To create strings

    GetResiduals();   // To calculate residual nucleus properties

    #ifdef debugQGSParticipants
      G4cout<<"All O.K. ======" <<G4endl;
    #endif
  }

  // clean-up, if necessary
  #ifdef debugQGSParticipants
     G4cout<<"Clearing "<<G4endl;
     G4cout<<"theInteractions.size() "<<theInteractions.size()<<G4endl;
  #endif

  if( Regge ) delete Regge; // Uzhi 18 Oct. 2016

  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();

// Erasing of target involved nucleons.
  #ifdef debugQGSParticipants
     G4cout<<"Erasing of involved target nucleons "<<NumberOfInvolvedNucleonsOfTarget<<G4endl;
  #endif

  if ( NumberOfInvolvedNucleonsOfTarget != 0 ) {
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
      if ( (aNucleon != 0 ) && (aNucleon->GetStatus() >= 1) ) delete aNucleon;
     }
  }

// Erasing of projectile involved nucleons.
  if ( NumberOfInvolvedNucleonsOfProjectile != 0 ) {
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
       G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfProjectile[i]->GetSplitableHadron();
       if ( aNucleon ) delete aNucleon;
     }
  }

  #ifdef debugQGSParticipants
     G4cout<<"Delition of target nucleons from soft interactions "<<theTargets.size()
           <<G4endl<<G4endl;
  #endif
  std::for_each(theTargets.begin(), theTargets.end(), DeleteSplitableHadron());
  theTargets.clear();

  if ( theProjectileSplitable ) {
    delete theProjectileSplitable;
    theProjectileSplitable = 0;
  }
}

//===========================================================
void G4QGSParticipants::PrepareInitialState( const G4ReactionProduct& thePrimary ) 
{
  // Clearing of the arrays
  // Erasing of the projectile
  G4InteractionContent* anIniteraction = theInteractions[0];
  G4VSplitableHadron* pProjectile = anIniteraction->GetProjectile();
  if( pProjectile ) delete pProjectile;

  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();

  // Erasing of the envolved nucleons and target nucleons from diffraction dissociations
  theNucleus->StartLoop();
  G4Nucleon* aNucleon;
  while ( ( aNucleon = theNucleus->GetNextNucleon() ) ) 
  {
   if ( aNucleon->AreYouHit() ) {
     G4VSplitableHadron* splaNucleon = aNucleon->GetSplitableHadron();
     if ( (splaNucleon != 0) && (splaNucleon->GetStatus() >=1) ) delete splaNucleon;
     aNucleon->Hit(nullptr);
     NumberOfInvolvedNucleonsOfTarget--;
   } 
  }

  // Erasing of nuclear nucleons participated in soft interactions
  std::for_each(theTargets.begin(), theTargets.end(), DeleteSplitableHadron());
  theTargets.clear();

  // Preparation to a new attempt
  theProjectile = thePrimary;

  theNucleus->Init(theNucleus->GetMassNumber(), theNucleus->GetCharge());
  theNucleus->SortNucleonsIncZ();       
  DoLorentzBoost(-theCurrentVelocity);  // Lorentz boost of the target nucleus

  G4LorentzVector Tmp( 0.0, 0.0, 0.0, 0.0 );
  NumberOfInvolvedNucleonsOfTarget= 0;
  TargetResidualMassNumber       = theNucleus->GetMassNumber();
  TargetResidualCharge           = theNucleus->GetCharge();
  TargetResidualExcitationEnergy = 0.0;

  G4Nucleon* NuclearNucleon;
  while ( ( NuclearNucleon = theNucleus->GetNextNucleon() ) )
              {Tmp+=NuclearNucleon->Get4Momentum();} 

  TargetResidual4Momentum        = Tmp;
}

//===========================================================
void G4QGSParticipants::GetList( const G4ReactionProduct& thePrimary ) { 
  #ifdef debugQGSParticipants
     G4cout<<G4endl<<"G4QGSParticipants::GetList +++++++++++++"<<G4endl;
  #endif

// Uzhi Direction: True - Proj, False - Target
  theProjectileSplitable = new G4QGSMSplitableHadron(thePrimary, TRUE);
  theProjectileSplitable->SetStatus(1);

  G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
  G4LorentzVector aNucleonMomentum(0.,0.,0., 938.0*MeV);

  G4double SS=(aPrimaryMomentum + aNucleonMomentum).mag2();

  Regge->SetS(SS); // Uzhi 18 Oct. 2016

//--------------------------------------
  theNucleus->StartLoop();
  G4Nucleon * tNucleon = theNucleus->GetNextNucleon();

  if ( ! tNucleon ) {
    #ifdef debugQGSParticipants
       G4cout << "QGSM - BAD situation: pNucleon is NULL ! Leaving immediately!" << G4endl;
    #endif
    return;
  }

// Determination of participating nucleons of nucleus ------------------------------------

  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();

  G4int totalCuts = 0;
  G4int MaxPower=thePrimary.GetMomentum().mag()/(3.3*GeV); if(MaxPower < 1) MaxPower=1; // Uzhi 3.3 GeV ??? tune

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = -1;
  while( (theInteractions.size() == 0) && ++loopCounter < maxNumberOfLoops )
  {
   InteractionMode = ALL; // Mode = ALL, WITHOUT_R, NON_DIFF           // Uzhi 18 Oct. 2016

   // choose random impact parameter of a collision
   std::pair<G4double, G4double> theImpactParameter;

   theImpactParameter = theNucleus->ChooseImpactXandY(theNucleus->GetOuterRadius()+theNucleonRadius);
   G4double impactX = theImpactParameter.first;
   G4double impactY = theImpactParameter.second;

   #ifdef debugQGSParticipants
     G4cout<<"InteractionMode "<<InteractionMode<<G4endl;
     G4cout<<"Impact parameter (fm ) "<<std::sqrt(sqr(impactX)+sqr(impactY))/fermi<<" "<<G4endl;
   #endif

   // loop over nucleons to find collisions
   theNucleus->StartLoop();
   G4int nucleonCount = -1;
   G4QGSParticipants_NPart = 0;

   G4double Power=MaxPower;    // Uzhi 2016

 // Uzhi 18 Oct. 2016 ModelMode = SOFT;//DIFFRACTIVE;
   while( (tNucleon = theNucleus->GetNextNucleon()) )
   {
    if(Power <= 0.) break; // #################################################################
    nucleonCount++;

    G4LorentzVector nucleonMomentum=tNucleon->Get4Momentum();

    G4double Distance2 = sqr(impactX - tNucleon->GetPosition().x()) +
                         sqr(impactY - tNucleon->GetPosition().y());

    G4double Pint(0.);                    // A probability of interaction at given impact parameter 
    G4double Pprd(0.), Ptrd(0.), Pdd(0.); // Probabilities of Proj. diffr., Target diffr., Double diffr. 
    G4double Pnd (0.), Pnvr(0.);          // Probabilities of non-diffr. and quark exchange  
    G4int    NcutPomerons(0);             // Number of cutted pomerons

    Regge->GetProbabilities(std::sqrt(Distance2), InteractionMode,
			      Pint, Pprd, Ptrd, Pdd, Pnd, Pnvr);
     #ifdef debugQGSParticipants
	G4cout<<"Nucleon & its impact parameter: "<<nucleonCount<<" "<<std::sqrt(Distance2)/fermi<<" (fm)"<<G4endl;
        G4cout<<"Probability of interaction:     "<<Pint<<G4endl;
	G4cout<<"Probability of PrD, TrD, DD:    "<<Pprd<<" "<<Ptrd<<" "<<Pdd<<G4endl;
	G4cout<<"Probability of NonDiff, QuarkExc.: "<<Pnd<<" "<<Pnvr<<" in inel. inter."<<G4endl;
     #endif

//Pint=1.; Pprd=0.; Ptrd=0.; Pdd=0.; Pnd=0.; Pnvr=1.;  InteractionMode=0;  // Uzhi for testing purposes
    if (Pint > G4UniformRand())
    {                             // An interaction is happend.

     G4double rndNumber = G4UniformRand();
     G4int InteractionType(0);

     if((InteractionMode==ALL)||(InteractionMode==WITHOUT_R))     // Mode = ALL, WITHOUT_R, NON_DIFF 
     {
	if( rndNumber < Pprd )                    {InteractionType = PrD;  InteractionMode = WITHOUT_R;}
	if( rndNumber < Pprd + Ptrd )             {InteractionType = TrD;  InteractionMode = WITHOUT_R;}
	if( rndNumber < Pprd + Ptrd + Pdd)        {InteractionType = DD;   InteractionMode = WITHOUT_R;}
	if( rndNumber < Pprd + Ptrd + Pdd + Pnd ) {InteractionType = NonD; InteractionMode = NON_DIFF;
						   NcutPomerons =  Regge->ncPomerons();                }
	if( InteractionMode == ALL ) 		  {InteractionType = Qexc;}
     }
     else  // InteractionMode == NON_DIFF
     {
	InteractionMode = NON_DIFF;
	if( rndNumber < Ptrd )                    {InteractionType = TrD; }
	if( rndNumber < Ptrd + Pnd)               {InteractionType = NonD;  NcutPomerons =  Regge->ncPomerons();}
     }

     if( (InteractionType == NonD) && (NcutPomerons == 0)) continue; 

     G4QGSParticipants_NPart ++;
     G4QGSMSplitableHadron* aTargetSPB = new G4QGSMSplitableHadron(*tNucleon);
     tNucleon->Hit(aTargetSPB);

     #ifdef debugQGSParticipants
	G4cout<<"An interaction is happend."<<G4endl;
        G4cout<<"Target nucleon - "<<nucleonCount<<" "
              <<tNucleon->GetDefinition()->GetParticleName()<<G4endl;
        G4cout<<"Interaction type:"<<InteractionType
              <<" (0 -PrD, 1 - TrD, 2 - DD, 3 - NonD, 4 - Qexc)"<<G4endl;
        G4cout<<"New Inter.  mode:"<<InteractionMode
              <<" (0 -ALL, 1 - WITHOUT_R, 2 - NON_DIFF)"<<G4endl;
        if( InteractionType == NonD )
          G4cout<<"Number of cutted pomerons: "<<NcutPomerons<<G4endl;
     #endif

     if((InteractionType == PrD) || (InteractionType == TrD) || (InteractionType == DD) ||
	(InteractionType == Qexc))
     {                                  // diffractive-like interaction occurs
      #ifdef debugQGSParticipants
         G4cout<<"Diffractive-like interaction occurs"<<G4endl;
      #endif

      G4InteractionContent * aInteraction = new G4InteractionContent(theProjectileSplitable);
      theProjectileSplitable->SetStatus(1*theProjectileSplitable->GetStatus());

      aInteraction->SetTarget(aTargetSPB);
      aInteraction->SetTargetNucleon(tNucleon);
      aTargetSPB->SetCollisionCount(0);
      aTargetSPB->SetStatus(1);

      aInteraction->SetNumberOfDiffractiveCollisions(1);
      aInteraction->SetNumberOfSoftCollisions(0);
      aInteraction->SetStatus(InteractionType);
      theInteractions.push_back(aInteraction);
     }
     else
     {                               // nondiffractive interaction occurs
        #ifdef debugQGSParticipants
           G4cout<<"Non-diffractive interaction occurs, max NcutPomerons "<<NcutPomerons<<G4endl;
        #endif

	G4int nCuts;

        G4int Vncut=0; 
	for(nCuts = 0; nCuts < NcutPomerons; nCuts++) 
	{       
	 if( G4UniformRand() < Power/MaxPower ){Vncut++; Power--; if(Power <= 0.) break;}
	}
        nCuts=Vncut;

//nCut=1;    // Uzhi for testing purposes
	if( nCuts == 0 ) {delete aTargetSPB; tNucleon->Hit(nullptr); continue;} 

        totalCuts += nCuts;                              // Uzhi 2016 ???
        #ifdef debugQGSParticipants
           G4cout<<"Number of cuts in the interaction "<<nCuts<<G4endl;
        #endif

	aTargetSPB->IncrementCollisionCount(nCuts);
        aTargetSPB->SetStatus(0);
        theTargets.push_back(aTargetSPB);

	theProjectileSplitable->IncrementCollisionCount(nCuts);
        theProjectileSplitable->SetStatus(0*theProjectileSplitable->GetStatus());

	G4InteractionContent * aInteraction = 
                                           new G4InteractionContent(theProjectileSplitable);
	aInteraction->SetTarget(aTargetSPB);
        aInteraction->SetTargetNucleon(tNucleon);
	aInteraction->SetNumberOfSoftCollisions(nCuts);
        aInteraction->SetStatus(InteractionType);
	theInteractions.push_back(aInteraction);
     }
    }    // End of if (Pint > G4UniformRand())
   }     // End of while( (tNucleon = theNucleus->GetNextNucleon()) )

  #ifdef debugQGSParticipants
    G4cout << G4endl<<"Number of wounded nucleons "<<G4QGSParticipants_NPart<<G4endl;
  #endif

  }  // End of while( (theInteractions.size() == 0) && ++loopCounter < maxNumberOfLoops )

  if ( loopCounter >= maxNumberOfLoops ) {
    #ifdef debugQGSParticipants
       G4cout <<"BAD situation: forced loop exit!" << G4endl;
    #endif
    // Perhaps there is something to set here...
    // Decrease impact parameter ??
    // Select collisions with only diffraction ??
    // Selecy only non-diffractive interactions ??
  }
//------------------------------------------------------------
  std::vector<G4InteractionContent*>::iterator i;

  if( InteractionMode == ALL )  // It can be if all interactions were quark-exchange. 
  {                             // Only the first one will be saved, all other will be erased.
     i = theInteractions.end()-1;

     while ( theInteractions.size() != 1 )  
     {
	G4InteractionContent* anInteraction = *i;
        G4Nucleon * pNucleon = anInteraction->GetTargetNucleon(); pNucleon->Hit(nullptr);

	delete *i;
	i=theInteractions.erase(i);
	i--;
     }          // End of while 
  }
  else
  {                             // All quark exchanges will be erased
     i = theInteractions.begin();
     while ( i != theInteractions.end() )  
     {
	G4InteractionContent* anInteraction = *i;

        if( anInteraction->GetStatus() == Qexc )
        {
         G4Nucleon*        aTargetNucleon = anInteraction->GetTargetNucleon();
	 aTargetNucleon->Hit(nullptr);

	 delete *i;
	 i=theInteractions.erase(i);
        }
        else
        {
         i++;
        }
     }          // End of while ( i != theInteractions.end() ) 
  }


  #ifdef debugQGSParticipants
    G4cout <<"Total number of cuts "<< totalCuts <<G4endl;
  #endif

}

//=============================================================
void G4QGSParticipants::StoreInvolvedNucleon() 
{ //To store nucleons involved in the interaction

  NumberOfInvolvedNucleonsOfTarget = 0;

  theNucleus->StartLoop();

  G4Nucleon* aNucleon;
  while ( ( aNucleon = theNucleus->GetNextNucleon() ) ) {
    if ( aNucleon->AreYouHit() ) {
      TheInvolvedNucleonsOfTarget[NumberOfInvolvedNucleonsOfTarget] = aNucleon;
      NumberOfInvolvedNucleonsOfTarget++;
    }
  }

  #ifdef debugQGSParticipants
     G4cout << G4endl<<"G4QGSParticipants::StoreInvolvedNucleon() if they were "<<G4endl
            <<"Stored # of wounded nucleons of target "
           << NumberOfInvolvedNucleonsOfTarget <<G4endl;
  #endif
  return;
}                        

//=============================================================

void G4QGSParticipants::ReggeonCascade() 
{ // Implementation of the reggeon theory inspired model of nuclear destruction 
  #ifdef debugQGSParticipants
     G4cout << G4endl<<"Reggeon cascading ........."<<G4endl;
     G4cout<<"C of nucl. desctruction "<<GetCofNuclearDestruction()
           <<" R2 "<<GetR2ofNuclearDestruction()/fermi/fermi<<" fermi^2"<<G4endl; 
  #endif

  G4int InitNINt = NumberOfInvolvedNucleonsOfTarget;

  // Reggeon cascading in target nucleus
  for ( G4int InvTN = 0; InvTN < InitNINt; InvTN++ ) { 
    G4Nucleon* aTargetNucleon = TheInvolvedNucleonsOfTarget[ InvTN ];

    G4double CreationTime = aTargetNucleon->GetSplitableHadron()->GetTimeOfCreation();

    G4double XofWoundedNucleon = aTargetNucleon->GetPosition().x();
    G4double YofWoundedNucleon = aTargetNucleon->GetPosition().y();
           
    G4V3DNucleus* theTargetNucleus = theNucleus;
    theTargetNucleus->StartLoop();

    G4int TrgNuc=0;
    G4Nucleon* Neighbour(0);
    while ( ( Neighbour = theTargetNucleus->GetNextNucleon() ) ) {
      TrgNuc++;
      if ( ! Neighbour->AreYouHit() ) {
        G4double impact2 = sqr( XofWoundedNucleon - Neighbour->GetPosition().x() ) +
                           sqr( YofWoundedNucleon - Neighbour->GetPosition().y() );

        if ( G4UniformRand() < GetCofNuclearDestruction() *
                               G4Exp( -impact2 / GetR2ofNuclearDestruction() )
           ) {  
          // The neighbour nucleon is involved in the reggeon cascade
          #ifdef debugQGSParticipants
             G4cout<<"Target nucleon involved in reggeon cascading No "<<TrgNuc<<" "
                   <<Neighbour->GetDefinition()->GetParticleName()<<G4endl;
          #endif
          TheInvolvedNucleonsOfTarget[ NumberOfInvolvedNucleonsOfTarget ] = Neighbour;
          NumberOfInvolvedNucleonsOfTarget++;

          G4QGSMSplitableHadron* targetSplitable = new G4QGSMSplitableHadron( *Neighbour ); 

          Neighbour->Hit( targetSplitable );
          targetSplitable->SetTimeOfCreation( CreationTime ); 
          targetSplitable->SetStatus( 2 );
          targetSplitable->SetCollisionCount(0);

          G4InteractionContent * anInteraction = new G4InteractionContent(theProjectileSplitable);
          anInteraction->SetTarget(targetSplitable);
          anInteraction->SetTargetNucleon(Neighbour);

          anInteraction->SetNumberOfDiffractiveCollisions(1);
          anInteraction->SetNumberOfSoftCollisions(0);
          anInteraction->SetStatus(3);                                                          // Uzhi (2); ???
          theInteractions.push_back(anInteraction);
        }
      } // End of if ( ! Neighbour->AreYouHit() )
    } // End of while ( ( Neighbour = theTargetNucleus->GetNextNucleon() ) )
  } // End of for ( G4int InvTN = 0; InvTN < InitNINt; InvTN++ )

  #ifdef debugQGSParticipants
     G4cout <<"Number of new involved nucleons "<<NumberOfInvolvedNucleonsOfTarget - InitNINt<<G4endl;
  #endif
  return;
}   

//============================================================================

G4bool G4QGSParticipants::PutOnMassShell() {

  G4bool isProjectileNucleus = false;
  if ( GetProjectileNucleus() ) {
    isProjectileNucleus = true;
  }

  #ifdef debugPutOnMassShell
     G4cout <<G4endl<< "PutOnMassShell start ..............." << G4endl;
     if ( isProjectileNucleus ) {G4cout << "PutOnMassShell for Nucleus_Nucleus " << G4endl;}
  #endif

  G4LorentzVector Pprojectile( theProjectile.GetMomentum(), theProjectile.GetTotalEnergy() );
  if ( Pprojectile.z() < 0.0 ) {
    return false;
  }

  G4bool isOk = true;
  
  G4LorentzVector Ptarget( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PtargetResidual( 0.0, 0.0, 0.0, 0.0 );
  G4double SumMasses = 0.0;
  G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
  G4double TargetResidualMass = 0.0; 

  #ifdef debugPutOnMassShell
     G4cout << "Target : ";
  #endif

  isOk = ComputeNucleusProperties( theTargetNucleus, Ptarget, PtargetResidual, SumMasses,
                                   TargetResidualExcitationEnergy, TargetResidualMass,
                                   TargetResidualMassNumber, TargetResidualCharge );
  if ( ! isOk ) return false;

  G4double Mprojectile  = 0.0;
  G4double M2projectile = 0.0;
  G4LorentzVector Pproj( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PprojResidual( 0.0, 0.0, 0.0, 0.0 );
  G4V3DNucleus* thePrNucleus = GetProjectileNucleus();
  G4double PrResidualMass = 0.0;

  if ( ! isProjectileNucleus ) {  // hadron-nucleus collision
    Mprojectile  = Pprojectile.mag();
    M2projectile = Pprojectile.mag2();
    SumMasses += Mprojectile + 20.0*MeV;
  } else {  // nucleus-nucleus or antinucleus-nucleus collision

    #ifdef debugPutOnMassShell
       G4cout << "Projectile : ";
    #endif

    isOk = ComputeNucleusProperties( thePrNucleus, Pproj, PprojResidual, SumMasses,
                                     ProjectileResidualExcitationEnergy, PrResidualMass,
                                     ProjectileResidualMassNumber, ProjectileResidualCharge );
    if ( ! isOk ) return false;
  }

  G4LorentzVector Psum = Pprojectile + Ptarget;   
  G4double SqrtS = Psum.mag();
  G4double     S = Psum.mag2();

  #ifdef debugPutOnMassShell
     G4cout << "Pproj "<<Pprojectile<<G4endl;
     G4cout << "Ptarg "<<Ptarget<<G4endl;
     G4cout << "Psum " << Psum/GeV << " GeV" << G4endl << "SqrtS " << SqrtS/GeV << " GeV" << G4endl
            << "SumMasses, PrResidualMass and TargetResidualMass " << SumMasses/GeV << " " 
            << PrResidualMass/GeV << " " << TargetResidualMass/GeV << " GeV" << G4endl;
     G4cout << "Ptar res. "<<PtargetResidual<<G4endl;
  #endif

  if ( SqrtS < SumMasses ) {
    return false;  // It is impossible to simulate after putting nuclear nucleons on mass-shell.
  }

  // Try to consider also the excitation energy of the residual nucleus, if this is
  // possible, with the available energy; otherwise, set the excitation energy to zero.

  G4double savedSumMasses = SumMasses;
  if ( isProjectileNucleus ) {
    SumMasses -= std::sqrt( sqr( PrResidualMass ) + PprojResidual.perp2() );
    SumMasses += std::sqrt( sqr( PrResidualMass + ProjectileResidualExcitationEnergy ) 
                            + PprojResidual.perp2() ); 
  }
  SumMasses -= std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
  SumMasses += std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy )
                          + PtargetResidual.perp2() );

  if ( SqrtS < SumMasses ) {
    SumMasses = savedSumMasses;
    if ( isProjectileNucleus ) {
      ProjectileResidualExcitationEnergy = 0.0;
    }
    TargetResidualExcitationEnergy = 0.0;
  }

  TargetResidualMass += TargetResidualExcitationEnergy;
  if ( isProjectileNucleus ) {
    PrResidualMass += ProjectileResidualExcitationEnergy;
  }

  #ifdef debugPutOnMassShell
     if ( isProjectileNucleus ) {
       G4cout << "PrResidualMass ProjResidualExcitationEnergy " << PrResidualMass/GeV << " "
	      << ProjectileResidualExcitationEnergy << " MeV" << G4endl;
     }
     G4cout << "TargetResidualMass TargetResidualExcitationEnergy " << TargetResidualMass/GeV << " GeV " 
            << TargetResidualExcitationEnergy << " MeV" << G4endl
            << "Sum masses " << SumMasses/GeV << G4endl;
  #endif

  // Sampling of nucleons what can transfer to delta-isobars
  if ( isProjectileNucleus  &&  thePrNucleus->GetMassNumber() != 1 ) {
      isOk = GenerateDeltaIsobar( SqrtS, NumberOfInvolvedNucleonsOfProjectile,
                                  TheInvolvedNucleonsOfProjectile, SumMasses );       
  }
  if ( theTargetNucleus->GetMassNumber() != 1 ) {
    isOk = isOk  &&
           GenerateDeltaIsobar( SqrtS, NumberOfInvolvedNucleonsOfTarget,
                                TheInvolvedNucleonsOfTarget, SumMasses );
  }
  if ( ! isOk ) return false;

  // Now we know that it is kinematically possible to produce a final state made
  // of the involved nucleons (or corresponding delta-isobars) and a residual nucleus.
  // We have to sample the kinematical variables which will allow to define the 4-momenta
  // of the final state. The sampled kinematical variables refer to the center-of-mass frame.
  // Notice that the sampling of the transverse momentum corresponds to take into account
  // Fermi motion.

  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp = toCms*Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) {  // "String" moving backwards in c.m.s., abort collision!
    return false; 
  }

  G4LorentzRotation toLab( toCms.inverse() );
  
  G4double YprojectileNucleus = 0.0;
  if ( isProjectileNucleus ) {
    Ptmp = toCms*Pproj;                      
    YprojectileNucleus = Ptmp.rapidity();
  }
  Ptmp = toCms*Ptarget;                      
  G4double YtargetNucleus = Ptmp.rapidity();

  // Ascribing of the involved nucleons Pt and Xminus
  G4double DcorP = 0.0;
  if ( isProjectileNucleus ) {
    DcorP = GetDofNuclearDestruction() / thePrNucleus->GetMassNumber();
  }
  G4double DcorT       = GetDofNuclearDestruction() / theTargetNucleus->GetMassNumber();
  G4double AveragePt2  = GetPt2ofNuclearDestruction();
  G4double maxPtSquare = GetMaxPt2ofNuclearDestruction();

  #ifdef debugPutOnMassShell
     if ( isProjectileNucleus ) {
       G4cout << "Y projectileNucleus " << YprojectileNucleus << G4endl;
     }
     G4cout << "Y targetNucleus     " << YtargetNucleus << G4endl 
            << "Dcor " << GetDofNuclearDestruction()
            << " DcorP DcorT " << DcorP << " " << DcorT << " AveragePt2 " << AveragePt2 << G4endl;
  #endif

  G4double M2proj = M2projectile;  // Initialization needed only for hadron-nucleus collisions
  G4double WplusProjectile = 0.0;
  G4double M2target = 0.0;
  G4double WminusTarget = 0.0;
  G4int NumberOfTries = 0;
  G4double ScaleFactor = 1.0;
  G4bool OuterSuccess = true;

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {  // while ( ! OuterSuccess )
    OuterSuccess = true;
    const G4int maxNumberOfTries = 1000;
    do {  // while ( SqrtS < Mprojectile + std::sqrt( M2target ) )
      NumberOfTries++;
      if ( NumberOfTries == 100*(NumberOfTries/100) ) {
        // After many tries, it is convenient to reduce the values of DcorP, DcorT and
        // AveragePt2, so that the sampled momenta (respectively, pz, and pt) of the
	// involved nucleons (or corresponding delta-isomers) are smaller, and therefore
        // it is more likely to satisfy the momentum conservation.
        ScaleFactor /= 2.0;
        DcorP       *= ScaleFactor;
        DcorT       *= ScaleFactor;
        AveragePt2  *= ScaleFactor;
      }
      if ( isProjectileNucleus ) {
        // Sampling of kinematical properties of projectile nucleons
        isOk = SamplingNucleonKinematics( AveragePt2, maxPtSquare, DcorP, 
                                          thePrNucleus, PprojResidual, 
                                          PrResidualMass, ProjectileResidualMassNumber,
                                          NumberOfInvolvedNucleonsOfProjectile, 
                                          TheInvolvedNucleonsOfProjectile, M2proj );
      }
      // Sampling of kinematical properties of target nucleons
      isOk = isOk  &&
             SamplingNucleonKinematics( AveragePt2, maxPtSquare, DcorT, 
                                        theTargetNucleus, PtargetResidual, 
                                        TargetResidualMass, TargetResidualMassNumber,
                                        NumberOfInvolvedNucleonsOfTarget, 
                                        TheInvolvedNucleonsOfTarget, M2target );

      #ifdef debugPutOnMassShell
      G4cout << "SqrtS, Mp+Mt, Mp, Mt " << SqrtS/GeV << " " 
             << ( std::sqrt( M2proj ) + std::sqrt( M2target) )/GeV << " "
             << std::sqrt( M2proj )/GeV << " " << std::sqrt( M2target )/GeV << G4endl;
      #endif

      if ( ! isOk ) return false;
    } while ( ( SqrtS < std::sqrt( M2proj ) + std::sqrt( M2target ) ) &&
              ++NumberOfTries < maxNumberOfTries );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( NumberOfTries >= maxNumberOfTries ) {
      return false;
    }
    if ( isProjectileNucleus ) {
      isOk = CheckKinematics( S, SqrtS, M2proj, M2target, YprojectileNucleus, true, 
                              NumberOfInvolvedNucleonsOfProjectile, 
                              TheInvolvedNucleonsOfProjectile,
                              WminusTarget, WplusProjectile, OuterSuccess );
    }
    isOk = isOk  &&
           CheckKinematics( S, SqrtS, M2proj, M2target, YtargetNucleus, false, 
                            NumberOfInvolvedNucleonsOfTarget, TheInvolvedNucleonsOfTarget,
                            WminusTarget, WplusProjectile, OuterSuccess );
    if ( ! isOk ) return false;
  } while ( ( ! OuterSuccess ) && 
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    return false;
  }

  // Now the sampling is completed, and we can determine the kinematics of the
  // whole system. This is done first in the center-of-mass frame, and then it is boosted
  // to the lab frame. The transverse momentum of the residual nucleus is determined as
  // the recoil of each hadron (nucleon or delta) which is emitted, i.e. in such a way
  // to conserve (by construction) the transverse momentum.

  if ( ! isProjectileNucleus ) {  // hadron-nucleus collision

    G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
    G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
    Pprojectile.setPz( Pzprojectile ); 
    Pprojectile.setE( Eprojectile );

    #ifdef debugPutOnMassShell
    G4cout << "Proj after in CMS " << Pprojectile/GeV <<" GeV"<< G4endl;
    #endif

    Pprojectile.transform( toLab );  
    theProjectile.SetMomentum( Pprojectile.vect() );
    theProjectile.SetTotalEnergy( Pprojectile.e() );

    if ( theProjectileSplitable ) theProjectileSplitable->Set4Momentum(Pprojectile);

    #ifdef debugPutOnMassShell
       G4cout << "Final proj. mom in Lab. " <<theProjectile.GetMomentum()/GeV<<" "
                                            <<theProjectile.GetTotalEnergy()/GeV<<" GeV"<<G4endl;
    #endif

  } else {  // nucleus-nucleus or antinucleus-nucleus collision

    isOk = FinalizeKinematics( WplusProjectile, true, toLab, PrResidualMass, 
                               ProjectileResidualMassNumber, NumberOfInvolvedNucleonsOfProjectile,
                               TheInvolvedNucleonsOfProjectile, ProjectileResidual4Momentum );

    #ifdef debugPutOnMassShell
    G4cout << "Projectile Residual4Momentum in CMS " << ProjectileResidual4Momentum/GeV <<" GeV"<< G4endl;
    #endif

    if ( ! isOk ) return false;

    ProjectileResidual4Momentum.transform( toLab );

    #ifdef debugPutOnMassShell
    G4cout << "Projectile Residual4Momentum in Lab " << ProjectileResidual4Momentum/GeV <<" GeV"<< G4endl;
    #endif

  }

  isOk = FinalizeKinematics( WminusTarget, false, toLab, TargetResidualMass, 
                             TargetResidualMassNumber, NumberOfInvolvedNucleonsOfTarget,
                             TheInvolvedNucleonsOfTarget, TargetResidual4Momentum );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in CMS " << TargetResidual4Momentum/GeV << " GeV "<< G4endl;
  #endif

  if ( ! isOk ) return false;

  TargetResidual4Momentum.transform( toLab );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in Lab " << TargetResidual4Momentum/GeV << " GeV "<< G4endl;
  #endif

  return true;

}

//============================================================================

G4ThreeVector G4QGSParticipants::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!

  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * 
                                        ( G4Exp( -maxPtSquare/AveragePt2 ) -1.0 ) );
  }
  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;

  return G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0.0 );    
}
//============================================================================

G4bool G4QGSParticipants::
ComputeNucleusProperties( G4V3DNucleus* nucleus,               // input parameter 
                          G4LorentzVector& nucleusMomentum,    // input & output parameter
                          G4LorentzVector& residualMomentum,   // input & output parameter
                          G4double& sumMasses,                 // input & output parameter
                          G4double& residualExcitationEnergy,  // input & output parameter
                          G4double& residualMass,              // input & output parameter
                          G4int& residualMassNumber,           // input & output parameter
                          G4int& residualCharge ) {            // input & output parameter

  // This method, which is called only by PutOnMassShell, computes some nucleus properties for:
  // -  either the target nucleus (which is never an antinucleus): this for any kind
  //    of hadronic interaction (hadron-nucleus, nucleus-nucleus, antinucleus-nucleus);
  // -  or the projectile nucleus or antinucleus: this only in the case of nucleus-nucleus
  //    or antinucleus-nucleus interaction.
  // This method assumes that the all the parameters have been initialized by the caller;
  // the action of this method consists in modifying all these parameters, except the
  // first one. The return value is "false" only in the case the pointer to the nucleus
  // is null.

  if ( ! nucleus ) return false;

  G4double ExcitationEPerWoundedNucleon = GetExcitationEnergyPerWoundedNucleon();

  // Loop over the nucleons of the nucleus. 
  // The nucleons that have been involved in the interaction (either from Glauber or
  // Reggeon Cascading) will be candidate to be emitted.
  // All the remaining nucleons will be the nucleons of the candidate residual nucleus.
  // The variable sumMasses is the amount of energy corresponding to:
  //     1. transverse mass of each involved nucleon
  //     2. 20.0*MeV separation energy for each involved nucleon
  //     3. transverse mass of the residual nucleus
  // In this first evaluation of sumMasses, the excitation energy of the residual nucleus
  // (residualExcitationEnergy, estimated by adding a constant value to each involved
  // nucleon) is not taken into account.
  G4Nucleon* aNucleon = 0;
  nucleus->StartLoop();
  while ( ( aNucleon = nucleus->GetNextNucleon() ) ) {  /* Loop checking, 07.08.2015, A.Ribon */
    nucleusMomentum += aNucleon->Get4Momentum();
    if ( aNucleon->AreYouHit() ) {  // Involved nucleons
      // Consider in sumMasses the nominal, i.e. on-shell, masses of the nucleons
      // (not the current masses, which could be different because the nucleons are off-shell).
      sumMasses += std::sqrt( sqr( aNucleon->GetDefinition()->GetPDGMass() ) 
                              +  aNucleon->Get4Momentum().perp2() );                     
      sumMasses += 20.0*MeV;  // Separation energy for a nucleon

//      residualExcitationEnergy += ExcitationEPerWoundedNucleon;                   // Uzhi April 2015
      residualExcitationEnergy += -ExcitationEPerWoundedNucleon*
                                   G4Log( G4UniformRand());                         // Uzhi April 2015
      residualMassNumber--;
      // The absolute value below is needed only in the case of anti-nucleus.
      residualCharge -= std::abs( G4int( aNucleon->GetDefinition()->GetPDGCharge() ) );
    } else {   // Spectator nucleons
      residualMomentum += aNucleon->Get4Momentum();
    }
  }
  #ifdef debugPutOnMassShell
  G4cout << "ExcitationEnergyPerWoundedNucleon " << ExcitationEPerWoundedNucleon <<" MeV"<<G4endl
         << "\t Residual Charge, MassNumber " << residualCharge << " " << residualMassNumber
         << G4endl << "\t Initial Momentum " << nucleusMomentum/GeV<<" GeV"
         << G4endl << "\t Residual Momentum   " << residualMomentum/GeV<<" GeV"<<G4endl;
  #endif
  residualMomentum.setPz( 0.0 ); 
  residualMomentum.setE( 0.0 );
  if ( residualMassNumber == 0 ) {
    residualMass = 0.0;
    residualExcitationEnergy = 0.0;
  } else {
    residualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
                     GetIonMass( residualCharge, residualMassNumber );
    if ( residualMassNumber == 1 ) {
      residualExcitationEnergy = 0.0;
    }
  }
  sumMasses += std::sqrt( sqr( residualMass ) + residualMomentum.perp2() );
  return true;
}


//============================================================================

G4bool G4QGSParticipants::
GenerateDeltaIsobar( const G4double sqrtS,                  // input parameter
                     const G4int numberOfInvolvedNucleons,  // input parameter
                     G4Nucleon* involvedNucleons[],         // input & output parameter
                     G4double& sumMasses ) {                // input & output parameter

  // This method, which is called only by PutOnMassShell, check whether is possible to
  // re-interpret some of the involved nucleons as delta-isobars:
  // - either by replacing a proton (2212) with a Delta+ (2214),
  // - or by replacing a neutron (2112) with a Delta0 (2114).
  // The on-shell mass of these delta-isobars is ~1232 MeV, so  ~292-294 MeV  heavier than
  // the corresponding nucleon on-shell mass. However  400.0*MeV  is considered to estimate
  // the max number of deltas compatible with the available energy.
  // The delta-isobars are considered with the same transverse momentum as their
  // corresponding nucleons.
  // This method assumes that all the parameters have been initialized by the caller;
  // the action of this method consists in modifying (eventually) involveNucleons and
  // sumMasses. The return value is "false" only in the case that the input parameters
  // have unphysical values.

  if ( sqrtS < 0.0  ||  numberOfInvolvedNucleons <= 0  ||  sumMasses < 0.0 ) return false;

  //const G4double ProbDeltaIsobar = 0.05;  // Uzhi 6.07.2012
  //const G4double ProbDeltaIsobar = 0.25;  // Uzhi 13.06.2013 
  const G4double probDeltaIsobar = 0.10;  // A.R. 07.08.2013

  G4int maxNumberOfDeltas = G4int( (sqrtS - sumMasses)/(400.0*MeV) );
  G4int numberOfDeltas = 0;

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    //G4cout << "i maxNumberOfDeltas probDeltaIsobar " << i << " " << maxNumberOfDeltas
    //       << " " << probDeltaIsobar << G4endl;
    if ( G4UniformRand() < probDeltaIsobar  &&  numberOfDeltas < maxNumberOfDeltas ) {
      numberOfDeltas++;
      if ( ! involvedNucleons[i] ) continue;
      G4VSplitableHadron* splitableHadron = involvedNucleons[i]->GetSplitableHadron();
      G4double massNuc = std::sqrt( sqr( splitableHadron->GetDefinition()->GetPDGMass() )
                                    + splitableHadron->Get4Momentum().perp2() );
      //AR The absolute value below is needed in the case of an antinucleus. 
      G4int pdgCode = std::abs( splitableHadron->GetDefinition()->GetPDGEncoding() );
      const G4ParticleDefinition* old_def = splitableHadron->GetDefinition();
      G4int newPdgCode = pdgCode/10; newPdgCode = newPdgCode*10 + 4; // Delta
      if ( splitableHadron->GetDefinition()->GetPDGEncoding() < 0 ) newPdgCode *= -1;
      const G4ParticleDefinition* ptr = 
        G4ParticleTable::GetParticleTable()->FindParticle( newPdgCode );
      splitableHadron->SetDefinition( ptr );
      G4double massDelta = std::sqrt( sqr( splitableHadron->GetDefinition()->GetPDGMass() )
                                      + splitableHadron->Get4Momentum().perp2() );
      //G4cout << i << " " << sqrtS/GeV << " " << sumMasses/GeV << " " << massDelta/GeV
      //       << " " << massNuc << G4endl;
      if ( sqrtS < sumMasses + massDelta - massNuc ) {  // Change cannot be accepted!
        splitableHadron->SetDefinition( old_def );
        break;
      } else {  // Change is accepted
        sumMasses += ( massDelta - massNuc );
      }
    } 
  }
  //G4cout << "maxNumberOfDeltas numberOfDeltas " << maxNumberOfDeltas << " " 
  //       << numberOfDeltas << G4endl;
  return true;
}


//============================================================================

G4bool G4QGSParticipants::
SamplingNucleonKinematics( G4double averagePt2,                   // input parameter
                           const G4double maxPt2,                 // input parameter
                           G4double dCor,                         // input parameter
                           G4V3DNucleus* nucleus,                 // input parameter
                           const G4LorentzVector& pResidual,      // input parameter
                           const G4double residualMass,           // input parameter
                           const G4int residualMassNumber,        // input parameter
                           const G4int numberOfInvolvedNucleons,  // input parameter 
                           G4Nucleon* involvedNucleons[],         // input & output parameter
                           G4double& mass2 ) {                    // output parameter

  // This method, which is called only by PutOnMassShell, does the sampling of:
  // -  either the target nucleons: this for any kind of hadronic interactions
  //    (hadron-nucleus, nucleus-nucleus, antinucleus-nucleus);
  // -  or the projectile nucleons or antinucleons: this only in the case of
  //    nucleus-nucleus or antinucleus-nucleus interactions, respectively.
  // This method assumes that all the parameters have been initialized by the caller;
  // the action of this method consists in changing the properties of the nucleons
  // whose pointers are in the vector involvedNucleons, as well as changing the
  // variable mass2.

  if ( ! nucleus ) return false;

  if ( residualMassNumber == 0  &&  numberOfInvolvedNucleons == 1 ) {
    dCor = 0.0; 
    averagePt2 = 0.0;
  } 

  G4bool success = true;                            

  G4double SumMasses = residualMass; 
  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    SumMasses += aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
  }

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {  // while ( ! success )

    success = true;
    G4ThreeVector ptSum( 0.0, 0.0, 0.0 );
    G4double xSum = 0.0;

    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4ThreeVector tmpPt = GaussianPt( averagePt2, maxPt2 );
      ptSum += tmpPt;
      G4ThreeVector tmpX = GaussianPt( dCor*dCor, 1.0 );
      G4double x = tmpX.x() +
                   aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()/SumMasses;
      if ( x < 0.0  ||  x > 1.0 ) { 
        success = false; 
        break;
      }
      xSum += x;
      //AR The energy is in the lab (instead of cms) frame but it will not be used.
      G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }

    if ( xSum < 0.0  ||  xSum > 1.0 ) success = false;

    if ( ! success ) continue;

    G4double deltaPx = ( ptSum.x() - pResidual.x() ) / numberOfInvolvedNucleons;
    G4double deltaPy = ( ptSum.y() - pResidual.y() ) / numberOfInvolvedNucleons;
    G4double delta = 0.0;
    if ( residualMassNumber == 0 ) {
      delta = ( xSum - 1.0 ) / numberOfInvolvedNucleons;
    } else {
      delta = 0.0;
    }

    xSum = 1.0;
    mass2 = 0.0;
    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4double x = aNucleon->Get4Momentum().pz() - delta;
      xSum -= x;               
      if ( residualMassNumber == 0 ) {
        if ( x <= 0.0  ||  x > 1.0 ) {
          success = false; 
          break;
        }
      } else {
        if ( x <= 0.0  ||  x > 1.0  ||  xSum <= 0.0  ||  xSum > 1.0 ) {
          success = false; 
          break;
        }
      }                                          
      G4double px = aNucleon->Get4Momentum().px() - deltaPx;
      G4double py = aNucleon->Get4Momentum().py() - deltaPy;
      mass2 += ( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() )
                    + sqr( px ) + sqr( py ) ) / x;
      G4LorentzVector tmp( px, py, x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }

    if ( success  &&  residualMassNumber != 0 ) {
      mass2 += ( sqr( residualMass ) + pResidual.perp2() ) / xSum;
    }

    #ifdef debugPutOnMassShell
    G4cout << "success " << success << G4endl << " Mt " << std::sqrt( mass2 )/GeV << G4endl;
    #endif

  } while ( ( ! success ) && 
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    success = false;
  }

  return success;
}


//============================================================================

G4bool G4QGSParticipants::
CheckKinematics( const G4double sValue,                 // input parameter
                 const G4double sqrtS,                  // input parameter
                 const G4double projectileMass2,        // input parameter
                 const G4double targetMass2,            // input parameter
                 const G4double nucleusY,               // input parameter
                 const G4bool isProjectileNucleus,      // input parameter
                 const G4int numberOfInvolvedNucleons,  // input parameter 
                 G4Nucleon* involvedNucleons[],         // input parameter
                 G4double& targetWminus,                // output parameter
                 G4double& projectileWplus,             // output parameter
                 G4bool& success ) {                    // input & output parameter

  // This method, which is called only by PutOnMassShell, checks whether the
  // kinematics is acceptable or not.
  // This method assumes that all the parameters have been initialized by the caller;
  // notice that the input boolean parameter isProjectileNucleus is meant to be true
  // only in the case of nucleus or antinucleus projectile.
  // The action of this method consists in computing targetWminus and projectileWplus
  // and setting the parameter success to false in the case that the kinematics should
  // be rejeted.

  G4double decayMomentum2 = sqr( sValue ) + sqr( projectileMass2 ) + sqr( targetMass2 )
                            - 2.0*sValue*projectileMass2 - 2.0*sValue*targetMass2 
                            - 2.0*projectileMass2*targetMass2;
  targetWminus = ( sValue - projectileMass2 + targetMass2 + std::sqrt( decayMomentum2 ) )
                 / 2.0 / sqrtS;
  projectileWplus = sqrtS - targetMass2/targetWminus;
  G4double projectilePz = projectileWplus/2.0 - projectileMass2/2.0/projectileWplus;
  G4double projectileE  = projectileWplus/2.0 + projectileMass2/2.0/projectileWplus;
  G4double projectileY(1.0e5);
  if(projectileE - projectilePz > 0.) {                                  // Uzhi 20.05.2015
           projectileY  = 0.5 * G4Log( (projectileE + projectilePz)/
                                       (projectileE - projectilePz) );
  }
  G4double targetPz = -targetWminus/2.0 + targetMass2/2.0/targetWminus;
  G4double targetE  =  targetWminus/2.0 + targetMass2/2.0/targetWminus;
  G4double targetY  = 0.5 * G4Log( (targetE + targetPz)/(targetE - targetPz) );

  #ifdef debugPutOnMassShell
  G4cout << "decayMomentum2 " << decayMomentum2 << G4endl 
         << "\t targetWminus projectileWplus " << targetWminus << " " << projectileWplus << G4endl
         << "\t projectileY targetY " << projectileY << " " << targetY << G4endl;
  #endif

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    G4double mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double x = tmp.z();
    G4double pz = -targetWminus*x/2.0 + mt2/(2.0*targetWminus*x);
    G4double e =   targetWminus*x/2.0 + mt2/(2.0*targetWminus*x);
    if ( isProjectileNucleus ) {
      pz = projectileWplus*x/2.0 - mt2/(2.0*projectileWplus*x);
      e =  projectileWplus*x/2.0 + mt2/(2.0*projectileWplus*x);
    }
    G4double nucleonY = 0.5 * G4Log( (e + pz)/(e - pz) ); 

    #ifdef debugPutOnMassShell
    G4cout << "i nY pY nY-AY AY " << i << " " << nucleonY << " " << projectileY <<G4endl;
    #endif

    if ( std::abs( nucleonY - nucleusY ) > 2  ||  
         ( isProjectileNucleus  &&  targetY > nucleonY )  ||
         ( ! isProjectileNucleus  &&  projectileY < nucleonY ) ) {
      success = false; 
      break;
    } 
  }
  return true;
}  

  
//============================================================================

G4bool G4QGSParticipants::
FinalizeKinematics( const G4double w,                            // input parameter
                    const G4bool isProjectileNucleus,            // input parameter
                    const G4LorentzRotation& boostFromCmsToLab,  // input parameter
                    const G4double residualMass,                 // input parameter
                    const G4int residualMassNumber,              // input parameter
                    const G4int numberOfInvolvedNucleons,        // input parameter 
                    G4Nucleon* involvedNucleons[],               // input & output parameter
	            G4LorentzVector& residual4Momentum ) {       // output parameter

  // This method, which is called only by PutOnMassShell, finalizes the kinematics:
  // this method is called when we are sure that the sampling of the kinematics is
  // acceptable.
  // This method assumes that all the parameters have been initialized by the caller;
  // notice that the input boolean parameter isProjectileNucleus is meant to be true
  // only in the case of nucleus or antinucleus projectile: this information is needed
  // because the sign of pz (in the center-of-mass frame) in this case is opposite
  // with respect to the case of a normal hadron projectile.
  // The action of this method consists in modifying the momenta of the nucleons
  // (in the lab frame) and computing the residual 4-momentum (in the center-of-mass
  // frame).

  G4ThreeVector residual3Momentum( 0.0, 0.0, 1.0 );

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    residual3Momentum -= tmp.vect();
    G4double mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double x = tmp.z();
    G4double pz = -w * x / 2.0  +  mt2 / ( 2.0 * w * x );
    G4double e  =  w * x / 2.0  +  mt2 / ( 2.0 * w * x );
    // Reverse the sign of pz in the case of nucleus or antinucleus projectile
    if ( isProjectileNucleus ) pz *= -1.0;
    tmp.setPz( pz ); 
    tmp.setE( e );
    tmp.transform( boostFromCmsToLab );
    aNucleon->SetMomentum( tmp );
    G4VSplitableHadron* splitableHadron = aNucleon->GetSplitableHadron();
    splitableHadron->Set4Momentum( tmp );
    #ifdef debugPutOnMassShell                                   // Uzhi 14.05.2015
       G4cout << "Target involved nucleon No, name, 4Mom " 
              << i<<" "<<aNucleon->GetDefinition()->GetParticleName()<<" "<<tmp<< G4endl;
    #endif
  }

  G4double residualMt2 = sqr( residualMass ) + sqr( residual3Momentum.x() )
                       + sqr( residual3Momentum.y() );

  #ifdef debugPutOnMassShell
  G4cout <<G4endl<< "w residual3Momentum.z() " << w << " " << residual3Momentum.z() << G4endl;
  #endif

  G4double residualPz = 0.0;
  G4double residualE  = 0.0;
  if ( residualMassNumber != 0 ) {
    residualPz = -w * residual3Momentum.z() / 2.0 + 
                  residualMt2 / ( 2.0 * w * residual3Momentum.z() );
    residualE  =  w * residual3Momentum.z() / 2.0 + 
                  residualMt2 / ( 2.0 * w * residual3Momentum.z() );
    // Reverse the sign of residualPz in the case of nucleus or antinucleus projectile
    if ( isProjectileNucleus ) residualPz *= -1.0;
  }

  residual4Momentum.setPx( residual3Momentum.x() );
  residual4Momentum.setPy( residual3Momentum.y() );
  residual4Momentum.setPz( residualPz ); 
  residual4Momentum.setE( residualE );

  return true;
}

//======================================================
void G4QGSParticipants::PerformDiffractiveCollisions()
{
	#ifdef debugQGSParticipants
	   G4cout<<G4endl<<"PerformDiffractiveCollisions()......"<<G4endl
                        <<"theInteractions.size() "<<theInteractions.size()<<G4endl;
	#endif

	unsigned int i;
	for(i = 0; i < theInteractions.size(); i++)
	{
		G4InteractionContent* anIniteraction = theInteractions[i];
		#ifdef debugQGSParticipants
		   G4cout<<"Interaction # and its status "
		         <<i<<" "<<theInteractions[i]->GetStatus()<<G4endl;
		#endif

		G4int InterStatus = theInteractions[i]->GetStatus(); 
                if( (InterStatus == PrD) || (InterStatus == TrD) || (InterStatus == DD))
                {		// Selection of diffractive interactions
		   #ifdef debugQGSParticipants
		     G4cout<<"Simulation of diffractive interaction. "<<InterStatus<<" PrD/TrD/DD/ND/Qech - 0,1,2,3,4"<<G4endl;
		   #endif

		   G4VSplitableHadron* aTarget = anIniteraction->GetTarget();

		   #ifdef debugQGSParticipants
			G4cout<<"The proj. before inter "
                              <<theProjectileSplitable->Get4Momentum()<<" "
                              <<theProjectileSplitable->Get4Momentum().mag()<<G4endl;
			G4cout<<"The targ. before inter "
                              <<aTarget->Get4Momentum()<<" "
                              <<aTarget->Get4Momentum().mag()<<G4endl;
		   #endif

		   if( InterStatus == PrD ) 
			theSingleDiffExcitation.ExciteParticipants(theProjectileSplitable, aTarget, TRUE); 

		   if( InterStatus == TrD ) 
			theSingleDiffExcitation.ExciteParticipants(theProjectileSplitable, aTarget, FALSE);

		   if( InterStatus == DD ) 
			theDiffExcitaton.ExciteParticipants(theProjectileSplitable, aTarget);

		   #ifdef debugQGSParticipants
			G4cout<<"The proj. after  inter "
                              <<theProjectileSplitable->Get4Momentum()<<" "
                              <<theProjectileSplitable->Get4Momentum().mag()<<G4endl;
			G4cout<<"The targ. after  inter "
                              <<aTarget->Get4Momentum()<<" "
                              <<aTarget->Get4Momentum().mag()<<G4endl;
		   #endif
                }

                if( InterStatus == Qexc )
                {		// Quark exchange process
		   #ifdef debugQGSParticipants
		     G4cout<<"Simulation of interaction with quark exchange."<<G4endl;
		   #endif
		   G4VSplitableHadron* aTarget = anIniteraction->GetTarget();

		   #ifdef debugQGSParticipants
			G4cout<<"The proj. before inter "
                              <<theProjectileSplitable->Get4Momentum()<<" "
                              <<theProjectileSplitable->Get4Momentum().mag()<<G4endl;
			G4cout<<"The targ. before inter "
                              <<aTarget->Get4Momentum()<<" "
                              <<aTarget->Get4Momentum().mag()<<G4endl;
		   #endif

		   theQuarkExchange.ExciteParticipants(theProjectileSplitable, aTarget);

		   #ifdef debugQGSParticipants
			G4cout<<"The proj. after  inter "
                              <<theProjectileSplitable->Get4Momentum()<<" "
                              <<theProjectileSplitable->Get4Momentum().mag()<<G4endl;
			G4cout<<"The targ. after  inter "
                              <<aTarget->Get4Momentum()<<" "
                              <<aTarget->Get4Momentum().mag()<<G4endl;
		   #endif
                }
	}
}

//======================================================
G4bool G4QGSParticipants::DeterminePartonMomenta()
{
  if( ! theProjectileSplitable ) return false;

  const G4double aHugeValue = 1.0e+10;

  #ifdef debugQGSParticipants
     G4cout<<G4endl<<"DeterminePartonMomenta()......"<<G4endl;
     G4cout<<"theProjectile status (0 -nondiffr, #0 diffr./reggeon):  "<<theProjectileSplitable->GetStatus()<<G4endl;
  #endif

  if(theProjectileSplitable->GetStatus() != 0) {return false;} // There were only diffractive interactions.

  G4LorentzVector Projectile4Momentum  = theProjectileSplitable->Get4Momentum();
  G4LorentzVector Psum = Projectile4Momentum;

G4double VqM_pr(0.), VaqM_pr(0.), VqM_tr(350.), VqqM_tr(700);                                         // Uzhi 21 Sept. 2016
if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) != 0) {VqM_pr=350*MeV; VaqM_pr=700*MeV;}// Uzhi 21 Sept. 2016

//VqM_pr=0.; VaqM_pr=0.; VqM_tr=0.; VqqM_tr=0.;   // Uzhi for testing purposes
  #ifdef debugQGSParticipants
     G4cout<<"Projectile 4 momentum "<<Psum<<G4endl
           <<"Target nucleon momenta at start"<<G4endl;
  #endif

  std::vector<G4VSplitableHadron*>::iterator i;
  G4int NuclNo=0;

  for(i = theTargets.begin(); i != theTargets.end(); i++ )
  {
   Psum += (*i)->Get4Momentum();
   #ifdef debugQGSParticipants
      G4cout<<"Nusleus nucleon # and its 4Mom. "<<NuclNo<<" "<<(*i)->Get4Momentum()<<G4endl;
   #endif
   NuclNo++;
  }

  G4LorentzRotation toCms( -1*Psum.boostVector() );

  G4LorentzVector Ptmp = toCms*Projectile4Momentum;

  toCms.rotateZ( -1*Ptmp.phi() );
  toCms.rotateY( -1*Ptmp.theta() );
  G4LorentzRotation toLab(toCms.inverse());
  Projectile4Momentum.transform( toCms );
//  Ptarget.transform( toCms );

  #ifdef debugQGSParticipants
     G4cout<<G4endl<<"In CMS---------------"<<G4endl;
     G4cout<<"Projectile 4 Mom "<<Projectile4Momentum<<G4endl;
  #endif

  NuclNo=0;
  G4LorentzVector Target4Momentum(0.,0.,0.,0.);
  for(i = theTargets.begin(); i != theTargets.end(); i++ )
  {
   G4LorentzVector tmp= (*i)->Get4Momentum();  tmp.transform( toCms );
   (*i)->Set4Momentum( tmp );
   #ifdef debugQGSParticipants
      G4cout<<"Target nucleon # and 4Mom "<<" "<<NuclNo<<" "<<(*i)->Get4Momentum()<<G4endl;
   #endif
   Target4Momentum += tmp;
   NuclNo++;
  }

  G4double S     = Psum.mag2();
  G4double SqrtS = std::sqrt(S);

  #ifdef debugQGSParticipants
     G4cout<<"Sum of target nucleons 4 momentum "<<Target4Momentum<<G4endl<<G4endl;
     G4cout<<"Target nucleons mom: px, py, z_1, m_i"<<G4endl;
  #endif

//G4double PplusProjectile = Projectile4Momentum.plus();
  G4double PminusTarget    = Target4Momentum.minus();
  NuclNo=0;
  
  for(i = theTargets.begin(); i != theTargets.end(); i++ )
  {
   G4LorentzVector tmp = (*i)->Get4Momentum(); // tmp.boost(bstToCM);

   //AR-19Jan2017 : the following line is causing a strange crash when Geant4
   //               is built in optimized mode.
   //               To fix it, I get the mass square instead of directly the
   //               mass from the Lorentz vector, and then I take care of the
   //               square root. If the mass square is negative, a JustWarning
   //               exception is thrown, and the mass is set to 0.
   //G4double Mass = tmp.mag();
   G4double Mass2 = tmp.mag2();
   G4double Mass = 0.0;
   if ( Mass2 < 0.0 ) {
     G4ExceptionDescription ed;
     ed << "Projectile " << theProjectile.GetDefinition()->GetParticleName()
        << "  4-momentum " << Psum << G4endl;
     ed << "LorentzVector tmp " << tmp << "  with Mass2 " << Mass2 << G4endl;
     G4Exception( "G4QGSParticipants::DeterminePartonMomenta(): 4-momentum with negative mass!",
                  "HAD_QGSPARTICIPANTS_001", JustWarning, ed );
   } else {
     Mass = sqrt( Mass2 );
   }

   tmp.setPz(tmp.minus()/PminusTarget);   tmp.setE(Mass);
   (*i)->Set4Momentum(tmp); 
   #ifdef debugQGSParticipants
      G4cout<<"Target nucleons # and mom: "<<NuclNo<<" "<<(*i)->Get4Momentum()<<G4endl;
   #endif
   NuclNo++;
  }

//+++++++++++++++++++++++++++++++++++++++++++
//G4double sigmaPt=0.5*GeV;                      // Uzhi 2016

  G4double SigPt = sigmaPt;   // Uzhi for testing purposes = 0
  G4Parton* aParton(0);
  G4ThreeVector aPtVector(0.,0.,0.);
  G4LorentzVector tmp(0.,0.,0.,0.);

  G4double Mt(0.);                             // Uzhi 19 Sept. 2016
  G4double ProjSumMt(0.), ProjSumMt2perX(0.);
  G4double TargSumMt(0.), TargSumMt2perX(0.);


  G4double aBeta = beta;   // Member of the class

  const G4ParticleDefinition* theProjectileDefinition = theProjectileSplitable->GetDefinition();
  if (theProjectileDefinition == G4PionMinus::PionMinusDefinition()) aBeta = 1.;  // Uzhi -0.5 
  if (theProjectileDefinition == G4Gamma::GammaDefinition())         aBeta = 1.;
  if (theProjectileDefinition == G4PionPlus::PionPlusDefinition())   aBeta = 1.;  // Uzhi -0.5
  if (theProjectileDefinition == G4PionZero::PionZeroDefinition())   aBeta = 1.;  // Uzhi -0.5
  if (theProjectileDefinition == G4KaonPlus::KaonPlusDefinition())   aBeta = 0.;
  if (theProjectileDefinition == G4KaonMinus::KaonMinusDefinition()) aBeta = 0.;

  G4double Xmin = 0.; // ==========================

  G4bool Success = true;  G4int attempt = 0;
  const G4int maxNumberOfAttempts = 1000;                  // Uzhi ############################
  do                         // while(!Success)
  {
   attempt++;  if( attempt ==  100*(attempt/100) ) {SigPt/=2.;}


   ProjSumMt=0.; ProjSumMt2perX=0.;
   TargSumMt=0.; TargSumMt2perX=0.;

   Success = true;
   G4int nSeaPair = theProjectileSplitable->GetSoftCollisionCount()-1;
   #ifdef debugQGSParticipants
      G4cout<<"attempt ------------------------ "<<attempt<<G4endl;
      G4cout<<"nSeaPair of proj "<<nSeaPair<<G4endl;
   #endif

   G4double SumPx = 0.; //theProjectileSplitable->Get4Momentum().px() * (-1.);
   G4double SumPy = 0.; //theProjectileSplitable->Get4Momentum().py() * (-1.);
   G4double SumZ = 0.;
   G4int               NumberOfUnsampledSeaQuarks = 2*nSeaPair;

   G4double Qmass=0.;
   for (G4int aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
   {
     aParton = theProjectileSplitable->GetNextParton();   // for quarks
     #ifdef debugQGSParticipants
        G4cout<<"Sea quarks: "<<aSeaPair<<" "<<aParton->GetDefinition()->GetParticleName();
     #endif
     aPtVector = GaussianPt(SigPt, aHugeValue);
     tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
     SumPx += aPtVector.x();   SumPy += aPtVector.y();
                     Mt = std::sqrt(aPtVector.mag2()+sqr(Qmass));       // Uzhi 19 Sept. 2016
                     ProjSumMt += Mt; 

//Sampling of Z fraction
//G4cout<<" NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<" ";
     tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ)); 
                     SumZ += tmp.z();

                     NumberOfUnsampledSeaQuarks--;
                     ProjSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);

     aParton = theProjectileSplitable->GetNextAntiParton();   // for anti-quarks
     #ifdef debugQGSParticipants
        G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
        G4cout<<"              "<<tmp<<" "<<SumZ<<G4endl;
     #endif
     aPtVector = GaussianPt(SigPt, aHugeValue);
     tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
     SumPx += aPtVector.x();   SumPy += aPtVector.y();
                     Mt = std::sqrt(aPtVector.mag2()+sqr(Qmass));       // Uzhi 19 Sept. 2016
                     ProjSumMt += Mt; 

//Sampling of Z fraction
//G4cout<<"NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<G4endl;
     tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ)); 
                     SumZ += tmp.z();

                     NumberOfUnsampledSeaQuarks--;
                     ProjSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);
     #ifdef debugQGSParticipants
        G4cout<<"              "<<tmp<<" "<<SumZ<<G4endl;
     #endif
   } 

//For valence quark                                                  // Uzhi 19 Sept. 2016
   aParton = theProjectileSplitable->GetNextParton();   // for quarks
   #ifdef debugQGSParticipants
      G4cout<<"Val quark of Pr"<<" "<<aParton->GetDefinition()->GetParticleName();
   #endif
   aPtVector = GaussianPt(SigPt, aHugeValue);
   tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
   SumPx += aPtVector.x();   SumPy += aPtVector.y();
                   Mt = std::sqrt(aPtVector.mag2()+sqr(VqM_pr));
                   ProjSumMt += Mt;

//Sampling of Z fraction
//G4cout<<" NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<" ";
   tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ)); 
                   SumZ += tmp.z();

                   ProjSumMt2perX +=sqr(Mt)/tmp.pz();
                   tmp.setE(sqr(Mt));
                   aParton->Set4Momentum(tmp);

// For valence di-quark
   aParton = theProjectileSplitable->GetNextAntiParton();
   #ifdef debugQGSParticipants
      G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
      G4cout<<"              "<<tmp<<" "<<SumZ<<" (z-fraction)"<<G4endl;
   #endif
   tmp.setPx(-SumPx);                     tmp.setPy(-SumPy);
                   Mt = std::sqrt(aPtVector.mag2()+sqr(VaqM_pr));
                         ProjSumMt += Mt;
                         tmp.setPz(1.-SumZ);

                   ProjSumMt2perX +=sqr(Mt)/tmp.pz();  // QQmass=750 MeV
                         tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);
   #ifdef debugQGSParticipants
      G4cout<<"              "<<tmp<<" "<<SumZ+(1.-SumZ)<<" (z-fraction)"<<G4endl;
   #endif
// End of work with the projectile

// Work with target nucleons 

   NuclNo=0;
   for(i = theTargets.begin(); i != theTargets.end(); i++ )
   {
    nSeaPair = (*i)->GetSoftCollisionCount()-1;
    #ifdef debugQGSParticipants
       G4cout<<"nSeaPair of target N "<<nSeaPair<<G4endl
             <<"Target nucleon 4Mom "<<(*i)->Get4Momentum()<<G4endl;
    #endif
//----------------------
    SumPx = (*i)->Get4Momentum().px() * (-1.);
    SumPy = (*i)->Get4Momentum().py() * (-1.);
    SumZ  = 0.                                ;   //**********

    G4double SumZw=0.;
    NumberOfUnsampledSeaQuarks = 2*nSeaPair;

    Qmass=0;	
    for (G4int aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
    {
     aParton = (*i)->GetNextParton();   // for quarks
     #ifdef debugQGSParticipants
        G4cout<<"Sea quarks: "<<aSeaPair<<" "<<aParton->GetDefinition()->GetParticleName();
     #endif
     aPtVector = GaussianPt(SigPt, aHugeValue);
     tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
     SumPx += aPtVector.x();   SumPy += aPtVector.y();
                     Mt=std::sqrt(aPtVector.mag2()+sqr(Qmass));
                     TargSumMt += Mt; 

//Sampling of Z fraction
//G4cout<<" NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<" ";
     tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ));
                     SumZ += tmp.z();
          tmp.setPz((*i)->Get4Momentum().pz()*tmp.pz());
		     SumZw+=tmp.pz();
                     NumberOfUnsampledSeaQuarks--;
                     TargSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);

     aParton = (*i)->GetNextAntiParton();   // for anti-quarks
     #ifdef debugQGSParticipants
        G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
        G4cout<<"              "<<tmp<<" "<<SumZw<<" "<<SumZ<<G4endl;
     #endif
     aPtVector = GaussianPt(SigPt, aHugeValue);
     tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
     SumPx += aPtVector.x();   SumPy += aPtVector.y();
                     Mt=std::sqrt(aPtVector.mag2()+sqr(Qmass));
                     TargSumMt += Mt; 

//Sampling of Z fraction
//G4cout<<" NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<G4endl;
     tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ)); 
                     SumZ += tmp.z();
     tmp.setPz((*i)->Get4Momentum().pz()*tmp.pz());
		     SumZw+=tmp.pz();
                     NumberOfUnsampledSeaQuarks--;
                     TargSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);
     #ifdef debugQGSParticipants
        G4cout<<"              "<<tmp<<" "<<SumZw<<" "<<SumZ<<G4endl;
     #endif
    } 

// Valence quark
    aParton = (*i)->GetNextParton();   // for quarks
    #ifdef debugQGSParticipants
       G4cout<<"Val quark of Tr"<<" "<<aParton->GetDefinition()->GetParticleName();
    #endif
    aPtVector = GaussianPt(SigPt, aHugeValue);
    tmp.setPx(aPtVector.x()); tmp.setPy(aPtVector.y());
    SumPx += aPtVector.x();   SumPy += aPtVector.y();
                     Mt=std::sqrt(aPtVector.mag2()+sqr(VqM_tr));
                     TargSumMt += Mt; 

//Sampling of Z fraction
//G4cout<<" NumberOfUnsampledSeaQuarks "<<NumberOfUnsampledSeaQuarks<<" ";
    tmp.setPz(SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta)*(1.0-SumZ)); 
                     SumZ += tmp.z();
                tmp.setPz((*i)->Get4Momentum().pz()*tmp.pz());
		     SumZw+=tmp.pz();
                     TargSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);

// Valence di-quark
    aParton = (*i)->GetNextAntiParton();   // for quarks
    #ifdef debugQGSParticipants
       G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
       G4cout<<"              "<<tmp<<" "<<SumZw<<" (sum z-fracs) "<<SumZ<<" (total z-sum) "<<G4endl;
    #endif
    tmp.setPx(-SumPx);                  tmp.setPy(-SumPy);
                     Mt=std::sqrt(aPtVector.mag2()+sqr(VqqM_tr));
                     TargSumMt += Mt; 

               tmp.setPz((*i)->Get4Momentum().pz()*(1.0 - SumZ));
		     SumZw+=tmp.pz();
                     TargSumMt2perX +=sqr(Mt)/tmp.pz();
                     tmp.setE(sqr(Mt));
                     aParton->Set4Momentum(tmp);
    #ifdef debugQGSParticipants
       G4cout<<"              "<<tmp<<" "<<SumZw<<" "<<1.0<<" "<<(*i)->Get4Momentum().pz()<<G4endl;
    #endif

   }   // End of for(i = theTargets.begin(); i != theTargets.end(); i++ )

   if( ProjSumMt      + TargSumMt      > SqrtS ) {
		Success = false; continue;}
   if( std::sqrt(ProjSumMt2perX) + std::sqrt(TargSumMt2perX) > SqrtS ) {
		Success = false; continue;}

  } while( (!Success) &&
           attempt < maxNumberOfAttempts );  /* Loop checking, 07.08.2015, A.Ribon */

  if ( attempt >= maxNumberOfAttempts ) {
    return false;
  }
//+++++++++++++++++++++++++++++++++++++++++++

  G4double DecayMomentum2 =    sqr(S) + sqr(ProjSumMt2perX) + sqr(TargSumMt2perX)
               - 2.0*S*ProjSumMt2perX - 2.0*S*TargSumMt2perX - 2.0*ProjSumMt2perX*TargSumMt2perX;

  G4double targetWminus=( S - ProjSumMt2perX + TargSumMt2perX + std::sqrt( DecayMomentum2 ))/2.0/SqrtS;
  G4double projectileWplus = SqrtS - TargSumMt2perX/targetWminus;
//--------------------------------------
  G4LorentzVector Tmp(0.,0.,0.,0.);
  G4double z(0.);

  G4int nSeaPair = theProjectileSplitable->GetSoftCollisionCount()-1;
   #ifdef debugQGSParticipants
      G4cout<<"Backward transformation ===================="<<G4endl;
      G4cout<<"nSeaPair of proj "<<nSeaPair<<G4endl;
   #endif

  for (G4int aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
  {
     aParton = theProjectileSplitable->GetNextParton();   // for quarks
     #ifdef debugQGSParticipants
        G4cout<<"Sea quarks: "<<aSeaPair<<" "<<aParton->GetDefinition()->GetParticleName();
     #endif
     Tmp =aParton->Get4Momentum(); z=Tmp.z();

     Tmp.setPz(projectileWplus*z/2.0 - Tmp.e()/(2.0*z*projectileWplus));
     Tmp.setE( projectileWplus*z/2.0 + Tmp.e()/(2.0*z*projectileWplus)); 
     Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);

     aParton = theProjectileSplitable->GetNextAntiParton();          // for anti-quarks
     #ifdef debugQGSParticipants
        G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
        G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<G4endl;
     #endif
     Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
     Tmp.setPz(projectileWplus*z/2.0 - Tmp.e()/(2.0*z*projectileWplus));
     Tmp.setE( projectileWplus*z/2.0 + Tmp.e()/(2.0*z*projectileWplus)); 
     Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);
     #ifdef debugQGSParticipants
        G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<G4endl;
     #endif
  } 

//For valence quark
   aParton = theProjectileSplitable->GetNextParton();   // for quarks
   #ifdef debugQGSParticipants
      G4cout<<"Val quark of Pr"<<" "<<aParton->GetDefinition()->GetParticleName();
   #endif
   Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
   Tmp.setPz(projectileWplus*z/2.0 - Tmp.e()/(2.0*z*projectileWplus));
   Tmp.setE( projectileWplus*z/2.0 + Tmp.e()/(2.0*z*projectileWplus)); 
   Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);

// For valence di-quark
   aParton = theProjectileSplitable->GetNextAntiParton();
   #ifdef debugQGSParticipants
      G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
      G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<" (mass)"<<G4endl;
   #endif
   Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
   Tmp.setPz(projectileWplus*z/2.0 - Tmp.e()/(2.0*z*projectileWplus));
   Tmp.setE( projectileWplus*z/2.0 + Tmp.e()/(2.0*z*projectileWplus)); 
   Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);

   #ifdef debugQGSParticipants
      G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<" (mass)"<<G4endl;
   #endif
// End of work with the projectile

// Work with target nucleons 
   NuclNo=0;
   for(i = theTargets.begin(); i != theTargets.end(); i++ )
   {
    nSeaPair = (*i)->GetSoftCollisionCount()-1;
    #ifdef debugQGSParticipants
       G4cout<<"nSeaPair of target and N# "<<nSeaPair<<" "<<NuclNo<<G4endl;
    #endif
    NuclNo++;
    for (G4int aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
    {
     aParton = (*i)->GetNextParton();   // for quarks
     #ifdef debugQGSParticipants
        G4cout<<"Sea quarks: "<<aSeaPair<<" "<<aParton->GetDefinition()->GetParticleName();
     #endif
     Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
     Tmp.setPz(-targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus));
     Tmp.setE(  targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus)); 
     Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);

     aParton = (*i)->GetNextAntiParton();   // for quarks
     #ifdef debugQGSParticipants
        G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
        G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<G4endl;
     #endif
     Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
     Tmp.setPz(-targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus));
     Tmp.setE(  targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus)); 
     Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);
     #ifdef debugQGSParticipants
        G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<G4endl;
     #endif
    } 

// Valence quark

    aParton = (*i)->GetNextParton();   // for quarks
    #ifdef debugQGSParticipants
       G4cout<<"Val quark of Tr"<<" "<<aParton->GetDefinition()->GetParticleName();
    #endif
    Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
    Tmp.setPz(-targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus));
    Tmp.setE(  targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus)); 
    Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);

// Valence di-quark
    aParton = (*i)->GetNextAntiParton();   // for quarks
    #ifdef debugQGSParticipants
       G4cout<<" "<<aParton->GetDefinition()->GetParticleName()<<G4endl;
       G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<" (mass)"<<G4endl;
    #endif
    Tmp =aParton->Get4Momentum(); z=Tmp.z(); 
    Tmp.setPz(-targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus));
    Tmp.setE(  targetWminus*z/2.0 + Tmp.e()/(2.0*z*targetWminus)); 
    Tmp.transform( toLab );

                     aParton->Set4Momentum(Tmp);
    #ifdef debugQGSParticipants
       G4cout<<"              "<<Tmp<<" "<<Tmp.mag()<<" (mass)"<<G4endl;
    #endif
NuclNo++;
   }   // End of for(i = theTargets.begin(); i != theTargets.end(); i++ )

//--------------------------------------
  return true;
}

//======================================================
G4double G4QGSParticipants::
SampleX(G4double anXmin, G4int nSea, G4int totalSea, G4double aBeta)
{
 G4double Xmin=anXmin; G4int Nsea=totalSea;   Xmin*=1.; Nsea++;// Must be erased Uzhi 12.05.2015 
   G4double Oalfa = 1./(alpha + 1.);
   G4double Obeta = 1./(aBeta + (alpha + 1.)*nSea + 1.); // ???????????
 
   G4double Ksi1, Ksi2, r1, r2, r12;
   const G4int maxNumberOfLoops = 1000;
   G4int loopCounter = 0;
   do
   {
     Ksi1 = G4UniformRand(); r1 = G4Pow::GetInstance()->powA(Ksi1,Oalfa);
     Ksi2 = G4UniformRand(); r2 = G4Pow::GetInstance()->powA(Ksi2,Obeta); 
     r12=r1+r2;
   } while( ( r12 > 1.) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
   if ( loopCounter >= maxNumberOfLoops ) {
     return 0.5;  // Just an acceptable value, without any physics consideration.
   }

   G4double result = r1/r12;
   return result;
} 

//======================================================
void G4QGSParticipants::CreateStrings()
{

   #ifdef debugQGSParticipants
   G4cout<<"CreateStrings() ..................."<<G4endl;
   #endif

   if ( ! theProjectileSplitable ) {
     #ifdef debugQGSParticipants
     G4cout<<"BAD situation: theProjectileSplitable is NULL ! Returning immediately"<<G4endl;
     #endif
     return;
   }

   #ifdef debugQGSParticipants
      G4cout<<"theProjectileSplitable->GetStatus() "<<theProjectileSplitable->GetStatus()<<G4endl;
      G4LorentzVector str4Mom;
   #endif

  if( theProjectileSplitable->GetStatus() == 1 ) // The projectile has participated only in diffr. inter.
  {
//?? unused     G4double      CreationTime = theProjectileSplitable->GetTimeOfCreation();
     G4ThreeVector Position     = theProjectileSplitable->GetPosition();

     G4PartonPair * aPair = new G4PartonPair(theProjectileSplitable->GetNextParton(), 
                                             theProjectileSplitable->GetNextAntiParton(),
			 	              G4PartonPair::DIFFRACTIVE, G4PartonPair::TARGET);
     #ifdef debugQGSParticipants
	G4cout << "Pr. Diffr. String: Qs 4mom X " <<G4endl;
	G4cout << "              " << aPair->GetParton1()->GetPDGcode()   << " "
				   << aPair->GetParton1()->Get4Momentum() << " "
				   << aPair->GetParton1()->GetX()         << " " << G4endl;
	G4cout << "              " << aPair->GetParton2()->GetPDGcode()   << " "
				   << aPair->GetParton2()->Get4Momentum() << " "
				   << aPair->GetParton2()->GetX()         << " " << G4endl;

	str4Mom += aPair->GetParton1()->Get4Momentum();
	str4Mom += aPair->GetParton2()->Get4Momentum();
     #endif

     thePartonPairs.push_back(aPair);
  }

  G4int N_EnvTarg = NumberOfInvolvedNucleonsOfTarget;

  for ( G4int i = 0; i < N_EnvTarg; i++ ) { 
    G4Nucleon* aTargetNucleon = TheInvolvedNucleonsOfTarget[ i ];

    G4VSplitableHadron* HitTargetNucleon = aTargetNucleon->GetSplitableHadron();

    #ifdef debugQGSParticipants
       G4cout<<"Involved Nucleon # and its status "<<i<<" "<<HitTargetNucleon->GetStatus()<<G4endl;
    #endif    
    if( HitTargetNucleon->GetStatus() >= 1) // Create diffractive string
    {
//unused??     G4double      CreationTime = HitTargetNucleon->GetTimeOfCreation();
     G4ThreeVector Position     = HitTargetNucleon->GetPosition();

     G4PartonPair * aPair = new G4PartonPair(HitTargetNucleon->GetNextParton(), 
                                             HitTargetNucleon->GetNextAntiParton(),
			 	           G4PartonPair::DIFFRACTIVE, G4PartonPair::TARGET);
     #ifdef debugQGSParticipants
	G4cout << "Tr. Diffr. String: Qs 4mom X " <<G4endl;
	G4cout << "Diffr. String " << aPair->GetParton1()->GetPDGcode()   << " "
				   << aPair->GetParton1()->Get4Momentum() << " "
				   << aPair->GetParton1()->GetX()         << " " << G4endl;
	G4cout << "              " << aPair->GetParton2()->GetPDGcode()   << " "
				   << aPair->GetParton2()->Get4Momentum() << " "
				   << aPair->GetParton2()->GetX()         << " " << G4endl;

	str4Mom += aPair->GetParton1()->Get4Momentum();
	str4Mom += aPair->GetParton2()->Get4Momentum();
     #endif

     thePartonPairs.push_back(aPair);
    }  // End of if( HitTargetNucleon->GetStatus() >= 1)
  }    // End of for ( G4int i = 0; i < N_EnvTarg; i++ )

//-----------------------------------------
  #ifdef debugQGSParticipants
     G4cout<<"Strings created in soft interactions"<<G4endl;
  #endif 
     std::vector<G4InteractionContent*>::iterator i;
     G4int IntNo=0;    
     i = theInteractions.begin();
     while ( i != theInteractions.end() )  /* Loop checking, 07.08.2015, A.Ribon */
     {
	G4InteractionContent* anIniteraction = *i;
	G4PartonPair * aPair = NULL;

        #ifdef debugQGSParticipants
           G4cout<<"An interaction # and soft coll. # "<<IntNo<<" "
                 <<anIniteraction->GetNumberOfSoftCollisions()<<G4endl;
        #endif 
        IntNo++;
	if (anIniteraction->GetNumberOfSoftCollisions())
	{
	  G4VSplitableHadron* pProjectile = anIniteraction->GetProjectile();
	  G4VSplitableHadron* pTarget     = anIniteraction->GetTarget();

	  for (G4int j = 0; j < anIniteraction->GetNumberOfSoftCollisions(); j++)
	  {
	    aPair = new G4PartonPair(pTarget->GetNextParton(), pProjectile->GetNextAntiParton(),
						G4PartonPair::SOFT, G4PartonPair::TARGET);
	    #ifdef debugQGSParticipants
		G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode()   << " "
		 		      << aPair->GetParton1()->Get4Momentum() << " "
				      <<aPair->GetParton1()->Get4Momentum().mag()<<G4endl;
				      //<< aPair->GetParton1()->GetX()         << " " << G4endl;
		G4cout << "         " << aPair->GetParton2()->GetPDGcode()   << " "
				      << aPair->GetParton2()->Get4Momentum() << " "
				      <<aPair->GetParton2()->Get4Momentum().mag()<<G4endl;
				      //<< aPair->GetParton2()->GetX()         << " " << G4endl;

		str4Mom += aPair->GetParton1()->Get4Momentum();
		str4Mom += aPair->GetParton2()->Get4Momentum();
	    #endif

	    thePartonPairs.push_back(aPair);

	    aPair = new G4PartonPair(pProjectile->GetNextParton(), pTarget->GetNextAntiParton(),
						G4PartonPair::SOFT, G4PartonPair::PROJECTILE);
	    #ifdef debugQGSParticipants
		G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode()   << " "
				      << aPair->GetParton1()->Get4Momentum() << " "
				      <<aPair->GetParton1()->Get4Momentum().mag()<<G4endl;
				      //<< aPair->GetParton1()->GetX()         << " " << G4endl;
		G4cout << "         " << aPair->GetParton2()->GetPDGcode()   << " "
				      << aPair->GetParton2()->Get4Momentum() << " "
				      <<aPair->GetParton2()->Get4Momentum().mag()<<G4endl;
				      //<< aPair->GetParton2()->GetX()         << " " << G4endl;
	    #endif
	    #ifdef debugQGSParticipants
		str4Mom += aPair->GetParton1()->Get4Momentum();
		str4Mom += aPair->GetParton2()->Get4Momentum();
	    #endif

	    thePartonPairs.push_back(aPair);
	  }  // End of for (G4int j = 0; j < anIniteraction->GetNumberOfSoftCollisions(); j++)

	   delete *i;
	   i=theInteractions.erase(i);    // i now points to the next interaction
	} else 
        {
          i++;
        }
     }          // End of while ( i != theInteractions.end() ) 
     #ifdef debugQGSParticipants
	G4cout << "Sum of strings 4 momenta " << str4Mom << G4endl<<G4endl;
     #endif
}

//============================================================================

void G4QGSParticipants::GetResiduals() {
  // This method is needed for the correct application of G4PrecompoundModelInterface

  #ifdef debugQGSParticipants
  G4cout << "GetResiduals(): GetProjectileNucleus()? "
         <<  GetProjectileNucleus() << G4endl;
  #endif

    #ifdef debugQGSParticipants
    G4cout << "NumberOfInvolvedNucleonsOfTarget "<< NumberOfInvolvedNucleonsOfTarget << G4endl;
    #endif

    G4double DeltaExcitationE = TargetResidualExcitationEnergy / 
                                G4double( NumberOfInvolvedNucleonsOfTarget );
    G4LorentzVector DeltaPResidualNucleus = TargetResidual4Momentum /
                                            G4double( NumberOfInvolvedNucleonsOfTarget );

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfTarget[i];

      #ifdef debugQGSParticipants
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      G4cout << i << " Hit? " << aNucleon->AreYouHit() << " " << targetSplitable << G4endl;
      if ( targetSplitable ) G4cout << i << "Status " << targetSplitable->GetStatus() << G4endl;
      #endif

      G4LorentzVector tmp = -DeltaPResidualNucleus;
      aNucleon->SetMomentum( tmp );
      aNucleon->SetBindingEnergy( DeltaExcitationE );
    }

//------------------------------------- 
    if( TargetResidualMassNumber != 0 )
    {
     G4ThreeVector bstToCM =TargetResidual4Momentum.findBoostToCM();

     G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
     G4LorentzVector residualMomentum(0.,0.,0.,0.);
     G4Nucleon* aNucleon = 0;
     theTargetNucleus->StartLoop();
     while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 07.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum(); tmp.boost(bstToCM);
         aNucleon->SetMomentum(tmp);
         residualMomentum +=tmp;
       }
     }

     residualMomentum/=TargetResidualMassNumber;

     G4double Mass = TargetResidual4Momentum.mag();
     G4double SumMasses=0.;
  
     aNucleon = 0;
     theTargetNucleus->StartLoop();
     while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 07.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum() - residualMomentum;
         G4double E=std::sqrt(tmp.vect().mag2()+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         tmp.setE(E);  aNucleon->SetMomentum(tmp);
         SumMasses+=E;
       }
     }

     G4double Chigh=Mass/SumMasses; G4double Clow=0; G4double C;
     const G4int maxNumberOfLoops = 1000;
     G4int loopCounter = 0;
     do
     {
      C=(Chigh+Clow)/2.;

      SumMasses=0.;
      aNucleon = 0;
      theTargetNucleus->StartLoop();
      while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 07.08.2015, A.Ribon */
        if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum();
         G4double E=std::sqrt(tmp.vect().mag2()*sqr(C)+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         SumMasses+=E;
        }
      }

      if(SumMasses > Mass) {Chigh=C;}
      else                 {Clow =C;}

     } while( (Chigh-Clow > 0.01) &&
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
     if ( loopCounter >= maxNumberOfLoops ) {
       #ifdef debugQGSParticipants
       G4cout <<"BAD situation: forced loop exit!" << G4endl;
       #endif
       // Perhaps there is something to set here...
     } else {

       aNucleon = 0;
       theTargetNucleus->StartLoop();
       while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 07.08.2015, A.Ribon */
         if ( !aNucleon->AreYouHit() ) { 
           G4LorentzVector tmp=aNucleon->Get4Momentum()*C;
           G4double E=std::sqrt(tmp.vect().mag2()+
                                sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
           tmp.setE(E); tmp.boost(-bstToCM);  
           aNucleon->SetMomentum(tmp);     
         }
       }
     }

    }   // End of if( TargetResidualMassNumber != 0 )
//-------------------------------------
   
    #ifdef debugQGSParticipants
    G4cout << "End GetResiduals -----------------" << G4endl;
    #endif

}


//======================================================
void G4QGSParticipants::PerformSoftCollisions()              // Uzhi It is not used
{
	std::vector<G4InteractionContent*>::iterator i;
	G4LorentzVector str4Mom;
	i = theInteractions.begin();
	while ( i != theInteractions.end() )  /* Loop checking, 07.08.2015, A.Ribon */
	{
		G4InteractionContent* anIniteraction = *i;
		G4PartonPair * aPair = NULL;
		if (anIniteraction->GetNumberOfSoftCollisions())
		{
			G4VSplitableHadron* pProjectile = anIniteraction->GetProjectile();
			G4VSplitableHadron* pTarget     = anIniteraction->GetTarget();
			for (G4int j = 0; j < anIniteraction->GetNumberOfSoftCollisions(); j++)
			{
				aPair = new G4PartonPair(pTarget->GetNextParton(), pProjectile->GetNextAntiParton(),
						G4PartonPair::SOFT, G4PartonPair::TARGET);
				#ifdef debugQGSParticipants
				G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " "
						<< aPair->GetParton1()->Get4Momentum() << " "
						<< aPair->GetParton1()->GetX() << " " << G4endl;
				G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " "
						<< aPair->GetParton2()->Get4Momentum() << " "
						<< aPair->GetParton2()->GetX() << " " << G4endl;
				#endif
				#ifdef debugQGSParticipants
				str4Mom += aPair->GetParton1()->Get4Momentum();
				str4Mom += aPair->GetParton2()->Get4Momentum();
				#endif
				thePartonPairs.push_back(aPair);
				aPair = new G4PartonPair(pProjectile->GetNextParton(), pTarget->GetNextAntiParton(),
						G4PartonPair::SOFT, G4PartonPair::PROJECTILE);
				#ifdef debugQGSParticipants
				G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " "
						<< aPair->GetParton1()->Get4Momentum() << " "
						<< aPair->GetParton1()->GetX() << " " << G4endl;
				G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " "
						<< aPair->GetParton2()->Get4Momentum() << " "
						<< aPair->GetParton2()->GetX() << " " << G4endl;
				#endif
				#ifdef debugQGSParticipants
				str4Mom += aPair->GetParton1()->Get4Momentum();
				str4Mom += aPair->GetParton2()->Get4Momentum();
				#endif
				thePartonPairs.push_back(aPair);
			}
			delete *i;
			i=theInteractions.erase(i);    // i now points to the next interaction
		} else {
			i++;
		}
	}
	#ifdef debugQGSParticipants
		G4cout << " string 4 mom " << str4Mom << G4endl;
	#endif
}


//************************************************
G4VSplitableHadron* G4QGSParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary) 
{
	// Check reaction threshold  - goes to CheckThreshold

	theProjectileSplitable = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE
        theProjectileSplitable->SetStatus(1);

	G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
	G4LorentzVector aTargetNMomentum(0.,0.,0.,938.);

	if((!(aPrimaryMomentum.e()>-1)) && (!(aPrimaryMomentum.e()<1)) )
	{
		throw G4HadronicException(__FILE__, __LINE__,
				"G4GammaParticipants::SelectInteractions: primary nan energy.");
	}
	G4double S = (aPrimaryMomentum + aTargetNMomentum).mag2();
	G4double ThresholdMass = thePrimary.GetMass() + 938.;
	ModelMode = SOFT;
	if (sqr(ThresholdMass + ThresholdParameter) > S)
	{
		ModelMode = DIFFRACTIVE;
	}
	if (sqr(ThresholdMass + QGSMThreshold) > S) // thus only diffractive in cascade!
	{
		ModelMode = DIFFRACTIVE;
	}

	// first find the collisions HPW
	std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
	theInteractions.clear();
	G4int totalCuts = 0;

	#ifdef debug_G4GammaParticipants
		G4double eK = thePrimary.GetKineticEnergy()/GeV;
		G4int nucleonCount = theNucleus->GetMassNumber();
	#endif

	G4int theCurrent = G4int(theNucleus->GetMassNumber()*G4UniformRand());
        G4int NucleonNo=0;

        theNucleus->StartLoop();
        G4Nucleon * pNucleon = 0;//theNucleus->GetNextNucleon();

        while( (pNucleon = theNucleus->GetNextNucleon()) )  /* Loop checking, 07.08.2015, A.Ribon */
        {if(NucleonNo == theCurrent) break; NucleonNo++;} 

        if ( pNucleon ) {
	  G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);

          pNucleon->Hit(aTarget);

	  if ( (0.06 > G4UniformRand() &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ) )
	  {
      		G4InteractionContent * aInteraction = new G4InteractionContent(theProjectileSplitable);
      		theProjectileSplitable->SetStatus(1*theProjectileSplitable->GetStatus());

      		aInteraction->SetTarget(aTarget);
      		aInteraction->SetTargetNucleon(pNucleon);
      		aTarget->SetCollisionCount(0);
      		aTarget->SetStatus(1);

      		aInteraction->SetNumberOfDiffractiveCollisions(1);
      		aInteraction->SetNumberOfSoftCollisions(0);
      		aInteraction->SetStatus(1);                     // Uzhi ???

      		theInteractions.push_back(aInteraction);
		totalCuts += 1;
	  }
	  else
	  {
		// nondiffractive soft interaction occurs
		aTarget->IncrementCollisionCount(1);
	        aTarget->SetStatus(0);
        	theTargets.push_back(aTarget);

		theProjectileSplitable->IncrementCollisionCount(1);
        	theProjectileSplitable->SetStatus(0*theProjectileSplitable->GetStatus());

		G4InteractionContent * aInteraction = 
                                           new G4InteractionContent(theProjectileSplitable);
		aInteraction->SetTarget(aTarget);
        	aInteraction->SetTargetNucleon(pNucleon);
		aInteraction->SetNumberOfSoftCollisions(1);
        	aInteraction->SetStatus(3);                    // Uzhi (0); ???
		theInteractions.push_back(aInteraction);
		totalCuts += 1;
	  }
        }
	return theProjectileSplitable;   //aProjectile;
}
