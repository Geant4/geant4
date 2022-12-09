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
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4QGSMFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4FragmentingString.hh"
#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4HadronicParameters.hh"
#include "G4Pow.hh"

//#define debug_QGSMfragmentation 

// Class G4QGSMFragmentation 
//****************************************************************************************
 
G4QGSMFragmentation::G4QGSMFragmentation() 
{
    SigmaQT = 0.45 * GeV;

    MassCut = 0.35*GeV; 

    SetStrangenessSuppression((1.0 - 0.16)/2.);

    // Check if charmed and bottom hadrons are enabled: if this is the case, then
    // set the non-zero probabilities for c-cbar and b-bbar creation from the vacuum,
    // else set them to 0.0. If these probabilities are/aren't zero then charmed or bottom
    // hadrons can't/can be created during the string fragmentation of ordinary
    // (i.e. not heavy) projectile hadron nuclear reactions.
    if ( G4HadronicParameters::Instance()->EnableBCParticles() ) {
      SetProbCCbar(0.0002);  // According to O.I. Piskunova Yad. Fiz. 56 (1993) 1094; tuned by Uzhi Oct. 2022
      SetProbBBbar(5.0e-5);  // According to O.I. Piskunova Yad. Fiz. 56 (1993) 1094
    } else {
      SetProbCCbar(0.0);
      SetProbBBbar(0.0);
    }
    
    SetDiquarkSuppression(0.32);
    SetDiquarkBreakProbability(0.7);

    SetMinMasses();

    arho = 0.5;    // alpha_rho0
    aphi = 0.0;    // alpha_fi
    aJPs =-2.2;    // alpha_J/Psi
    aUps =-8.0;    // alpha_Y      ??? O. Piskunova Yad. Phys. 56 (1993) 1094.

    aksi =-1.0; 
    alft = 0.5;    // 2 * alpha'_R *<Pt^2>

    an    = -0.5 ; 
    ala   = -0.75; // an - arho/2 + aphi/2 
    alaC  =  an - arho/2.0 + aJPs/2.0;  
    alaB  =  an - arho/2.0 + aUps/2.0;  
    aXi   =  0.0;  // ??
    aXiC  =  0.0;  // ??
    aXiB  =  0.0;  // ??
    aXiCC =  0.0;  // ??
    aXiCB =  0.0;  // ??
    aXiBB =  0.0;  // ??

    SetFFq2q();
    SetFFq2qq();
    SetFFqq2q();
    SetFFqq2qq();
                         // d  u   s   c   b
    G4int Index[5][5] = { { 0, 1,  2,  3,  4 },    // d
                          { 1, 5,  6,  7,  8 },    // u
                          { 2, 6,  9, 10, 11 },    // s
                          { 3, 7, 10, 12, 13 },    // c
                          { 4, 8, 11, 13, 14 } };  // b
    for (G4int i = 0; i < 5; i++ ) {
      for ( G4int j = 0; j < 5; j++ ) { 
        IndexDiQ[i][j] = Index[i][j];
      } 
    };
}

G4QGSMFragmentation::~G4QGSMFragmentation()
{}

//----------------------------------------------------------------------------------------------------------

G4KineticTrackVector* G4QGSMFragmentation::FragmentString(const G4ExcitedString& theString)
{

  G4FragmentingString  aString(theString);
  SetMinimalStringMass(&aString);

  #ifdef debug_QGSMfragmentation
  G4cout<<G4endl<<"QGSM StringFragm: String Mass "
                             <<theString.Get4Momentum().mag()<<" Pz "
                             <<theString.Get4Momentum().pz()
                             <<"------------------------------------"<<G4endl;
  G4cout<<"String ends Direct "<<theString.GetLeftParton()->GetPDGcode()<<" "
                               <<theString.GetRightParton()->GetPDGcode()<<" "
                               <<theString.GetDirection()<< G4endl;
  G4cout<<"Left  mom "<<theString.GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"Right mom "<<theString.GetRightParton()->Get4Momentum()<<G4endl;
  G4cout<<"Check for Fragmentation "<<G4endl;
  #endif

  // Can no longer modify Parameters for Fragmentation.
  PastInitPhase=true;
	
  // Check if string has enough mass to fragment...
  G4KineticTrackVector * LeftVector=NULL;

  if ( !IsItFragmentable(&aString) ) {
     LeftVector=ProduceOneHadron(&theString);

     #ifdef debug_QGSMfragmentation
     if ( LeftVector != 0 ) G4cout<<"Non fragmentable - the string is converted to one hadron "<<G4endl;
     #endif

     if ( LeftVector == nullptr ) LeftVector = new G4KineticTrackVector;
     return LeftVector;
  }

  #ifdef debug_QGSMfragmentation
  G4cout<<"The string will be fragmented. "<<G4endl;
  #endif
	
  LeftVector = new G4KineticTrackVector;
  G4KineticTrackVector * RightVector=new G4KineticTrackVector;

  G4ExcitedString *theStringInCMS=CopyExcited(theString);
  G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();

  G4bool success=false, inner_sucess=true;
  G4int attempt=0;
  while ( !success && attempt++ < StringLoopInterrupt )  /* Loop checking, 07.08.2015, A.Ribon */
  {
                #ifdef debug_QGSMfragmentation
                G4cout<<"Loop_toFrag "<<theStringInCMS->GetLeftParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetRightParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetDirection()<< G4endl;
                #endif

		G4FragmentingString *currentString=new G4FragmentingString(*theStringInCMS);

		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		RightVector->clear();
		
		inner_sucess=true;  // set false on failure..
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = -1;
		while (! StopFragmenting(currentString) && ++loopCounter < maxNumberOfLoops )   /* Loop checking, 07.08.2015, A.Ribon */
		{  // Split current string into hadron + new string

                        #ifdef debug_QGSMfragmentation
                        G4cout<<"The string can fragment. "<<G4endl;;
                        #endif
			G4FragmentingString *newString=0;  // used as output from SplitUp...
			G4KineticTrack * Hadron=Splitup(currentString,newString);

			if ( Hadron != 0 ) 
			{
                           #ifdef debug_QGSMfragmentation
                           G4cout<<"Hadron prod at fragm. "<<Hadron->GetDefinition()->GetParticleName()<<G4endl;
                           #endif
                           // To close the production of hadrons at fragmentation stage
			   if ( currentString->GetDecayDirection() > 0 )
				   LeftVector->push_back(Hadron);
       			   else
	  			   RightVector->push_back(Hadron);

			   delete currentString;
			   currentString=newString;

			} else {

                           #ifdef debug_QGSMfragmentation
                           G4cout<<"abandon ... start from the beginning ---------------"<<G4endl;
                           #endif

			   // Abandon ... start from the beginning
			   if (newString) delete newString;
			   inner_sucess=false;
			   break;
			}
		}
                if ( loopCounter >= maxNumberOfLoops ) {
                  inner_sucess=false;
                }

		// Split current string into 2 final Hadrons
                #ifdef debug_QGSMfragmentation
                if( inner_sucess ) {
                  G4cout<<"Split remaining string into 2 final hadrons."<<G4endl;
                } else {
		  G4cout<<" New attempt to fragment string"<<G4endl;
		}
                #endif
                // To the close production of hadrons at last string decay
		if ( inner_sucess && 
		     SplitLast(currentString,LeftVector, RightVector) ) 
		{
			success=true;
		}
		delete currentString;
	}
	
	delete theStringInCMS;
	
	if ( ! success )
	{
		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		delete RightVector;
		return LeftVector;
	}
		
	// Join Left- and RightVector into LeftVector in correct order.
	while(!RightVector->empty())  /* Loop checking, 07.08.2015, A.Ribon */
	{
	    LeftVector->push_back(RightVector->back());
	    RightVector->erase(RightVector->end()-1);
	}
	delete RightVector;

	CalculateHadronTimePosition(theString.Get4Momentum().mag(), LeftVector);

	G4LorentzRotation toObserverFrame(toCms.inverse());

	for (size_t C1 = 0; C1 < LeftVector->size(); C1++)
	{
	   G4KineticTrack* Hadron = LeftVector->operator[](C1);
	   G4LorentzVector Momentum = Hadron->Get4Momentum();
	   Momentum = toObserverFrame*Momentum;
	   Hadron->Set4Momentum(Momentum);
	   G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
	   Momentum = toObserverFrame*Coordinate;
	   Hadron->SetFormationTime(Momentum.e());
	   G4ThreeVector aPosition(Momentum.vect());
	   Hadron->SetPosition(theString.GetPosition()+aPosition);
	}
	return LeftVector;
}

//----------------------------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::IsItFragmentable(const G4FragmentingString * string)
{
	return sqr( MinimalStringMass + MassCut ) < string->Mass2();
}

//----------------------------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::StopFragmenting(const G4FragmentingString * string)
{
	SetMinimalStringMass(string);
        if ( MinimalStringMass < 0.0 ) return true;

        G4double smass = string->Mass();
	G4double x = (string->IsAFourQuarkString()) ? 0.005*(smass - MinimalStringMass)
	  : 0.66e-6*(smass - MinimalStringMass)*(smass + MinimalStringMass);

        G4bool res = true;
        if(x > 0.0) {
          res = (x < 200.) ? (G4UniformRand() < G4Exp(-x)) : false;
	}
	return res;
}

//-----------------------------------------------------------------------------

G4KineticTrack * G4QGSMFragmentation::Splitup( G4FragmentingString *string, 
			                       G4FragmentingString *&newString )
{
       #ifdef debug_QGSMfragmentation
       G4cout<<G4endl;
       G4cout<<"Start SplitUP (G4VLongitudinalStringDecay) ========================="<<G4endl;
       G4cout<<"String partons: " <<string->GetLeftParton()->GetPDGEncoding()<<" "
                                  <<string->GetRightParton()->GetPDGEncoding()<<" "
             <<"Direction "       <<string->GetDecayDirection()<<G4endl;
       #endif

       //... random choice of string end to use for creating the hadron (decay)   
       G4int SideOfDecay = (G4UniformRand() < 0.5)? 1: -1;
       if (SideOfDecay < 0)
       {
	  string->SetLeftPartonStable();
       } else
       {
          string->SetRightPartonStable();
       }

       G4ParticleDefinition *newStringEnd;
       G4ParticleDefinition * HadronDefinition;
       if (string->DecayIsQuark())
       {
	  G4double ProbDqADq = GetDiquarkSuppress();

	  G4int NumberOfpossibleBaryons = 2;

	  if (string->GetLeftParton()->GetParticleSubType()  != "quark") NumberOfpossibleBaryons++; 
	  if (string->GetRightParton()->GetParticleSubType() != "quark") NumberOfpossibleBaryons++;

	  G4double ActualProb  = ProbDqADq ;
          ActualProb *= (1.0-G4Exp(2.0*(1.0 - string->Mass()/(NumberOfpossibleBaryons*1400.0))));

	  SetDiquarkSuppression(ActualProb); 

       	  HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);

	  SetDiquarkSuppression(ProbDqADq);
       } else {
          HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
       }      

       if ( HadronDefinition == NULL ) return NULL;

       #ifdef debug_QGSMfragmentation
       G4cout<<"The parton "<<string->GetDecayParton()->GetPDGEncoding()<<" "
             <<" produces hadron "<<HadronDefinition->GetParticleName()
             <<" and is transformed to "<<newStringEnd->GetPDGEncoding()<<G4endl;
       G4cout<<"The side of the string decay Left/Right (1/-1) "<<SideOfDecay<<G4endl;
       #endif
       // create new String from old, ie. keep Left and Right order, but replace decay

       newString=new G4FragmentingString(*string,newStringEnd); // To store possible
                                                                // quark containt of new string

       #ifdef debug_QGSMfragmentation
       G4cout<<"An attempt to determine its energy (SplitEandP)"<<G4endl;
       #endif
       G4LorentzVector* HadronMomentum=SplitEandP(HadronDefinition, string, newString);

       delete newString; newString=0;
	
       G4KineticTrack * Hadron =0;
       if ( HadronMomentum != 0 ) {

           #ifdef debug_QGSMfragmentation                     
           G4cout<<"The attempt was successful"<<G4endl;
           #endif
	   G4ThreeVector   Pos;
	   Hadron = new G4KineticTrack(HadronDefinition, 0,Pos, *HadronMomentum);

     	   newString=new G4FragmentingString(*string,newStringEnd,HadronMomentum);

	   delete HadronMomentum;
       }
       else
       {
         #ifdef debug_QGSMfragmentation
         G4cout<<"The attempt was not successful !!!"<<G4endl;
         #endif
       }

       #ifdef debug_VStringDecay
       G4cout<<"End SplitUP (G4VLongitudinalStringDecay) ====================="<<G4endl;
       #endif

       return Hadron;
}

//-----------------------------------------------------------------------------

G4ParticleDefinition *G4QGSMFragmentation::DiQuarkSplitup( G4ParticleDefinition* decay,
                                                           G4ParticleDefinition *&created )
{
   //... can Diquark break or not?
   if (G4UniformRand() < DiquarkBreakProb )  //... Diquark break
   {
      G4int stableQuarkEncoding = decay->GetPDGEncoding()/1000;
      G4int decayQuarkEncoding = (decay->GetPDGEncoding()/100)%10;

      if (G4UniformRand() < 0.5)
      {
         G4int Swap = stableQuarkEncoding;
         stableQuarkEncoding = decayQuarkEncoding;
         decayQuarkEncoding = Swap;
      }

      G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1;  // if we have a quark, we need antiquark

      G4double StrSup=GetStrangeSuppress();
      SetStrangenessSuppression((1.0 - 0.07)/2.);  // Prob qq->K qq' 0.07
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      SetStrangenessSuppression(StrSup);

      //... Build new Diquark
      G4int QuarkEncoding=QuarkPair.second->GetPDGEncoding();
      G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
      G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
      
      created = FindParticle(NewDecayEncoding);
      G4ParticleDefinition * decayQuark=FindParticle(decayQuarkEncoding);
      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decayQuark);

      DecayQuark = decay->GetPDGEncoding();
      NewQuark   = NewDecayEncoding;

      return had;

   } else {  //... Diquark does not break

      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1;  // if we have a diquark, we need quark)
      G4double StrSup=GetStrangeSuppress();  // for changing s-sbar production
      SetStrangenessSuppression((1.0 - 0.07)/2.);
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      SetStrangenessSuppression(StrSup);

      created = QuarkPair.second;

      DecayQuark = decay->GetPDGEncoding();
      NewQuark   = created->GetPDGEncoding();

      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
      return had;
   }
}

//-----------------------------------------------------------------------------------------

G4LorentzVector * G4QGSMFragmentation::SplitEandP(G4ParticleDefinition * pHadron,
                                                  G4FragmentingString * string,  
                                                  G4FragmentingString * NewString)
{
       G4double HadronMass = pHadron->GetPDGMass();

       SetMinimalStringMass(NewString);

       if ( MinimalStringMass < 0.0 ) return nullptr;

       #ifdef debug_QGSMfragmentation
       G4cout<<"G4QGSMFragmentation::SplitEandP "<<pHadron->GetParticleName()<<G4endl;
       G4cout<<"String 4 mom, String M "<<string->Get4Momentum()<<" "<<string->Mass()<<G4endl;
       G4cout<<"HadM MinimalStringMassLeft StringM hM+sM "<<HadronMass<<" "<<MinimalStringMass<<" "
             <<string->Mass()<<" "<<HadronMass+MinimalStringMass<<G4endl;
       #endif

       if (HadronMass + MinimalStringMass > string->Mass())
       {
         #ifdef debug_QGSMfragmentation
         G4cout<<"Mass of the string is not sufficient to produce the hadron!"<<G4endl;
         #endif
	 return 0;
       }  // have to start all over!

       // calculate and assign hadron transverse momentum component HadronPx andHadronPy
       G4double StringMT2 = string->MassT2();
       G4double StringMT  = std::sqrt(StringMT2);

       G4LorentzVector String4Momentum = string->Get4Momentum();
       String4Momentum.setPz(0.);
       G4ThreeVector StringPt = String4Momentum.vect();

       G4ThreeVector HadronPt    , RemSysPt;
       G4double      HadronMassT2, ResidualMassT2;

       // Mt distribution is implemented
       G4double HadronMt, Pt, Pt2, phi;

       //...  sample Pt of the hadron
       G4int attempt=0;
       do
       {
         attempt++; if (attempt > StringLoopInterrupt) return 0;

         HadronMt = HadronMass - 200.0*G4Log(G4UniformRand());    // 200.0 must be tuned
         Pt2 = sqr(HadronMt)-sqr(HadronMass); Pt=std::sqrt(Pt2);
         phi = 2.*pi*G4UniformRand();
         G4ThreeVector SampleQuarkPtw= G4ThreeVector(Pt*std::cos(phi), Pt*std::sin(phi), 0);
         HadronPt =SampleQuarkPtw  + string->DecayPt();
         HadronPt.setZ(0);
         RemSysPt = StringPt - HadronPt;

         HadronMassT2 = sqr(HadronMass) + HadronPt.mag2();
         ResidualMassT2=sqr(MinimalStringMass) + RemSysPt.mag2();

       } while (std::sqrt(HadronMassT2) + std::sqrt(ResidualMassT2) > StringMT);  /* Loop checking, 07.08.2015, A.Ribon */

       //...  sample z to define hadron longitudinal momentum and energy
       //... but first check the available phase space

       G4double Pz2 = (sqr(StringMT2 - HadronMassT2 - ResidualMassT2) -
		      4*HadronMassT2 * ResidualMassT2)/4./StringMT2;

       if ( Pz2 < 0 ) {return 0;}          // have to start all over!

       //... then compute allowed z region  z_min <= z <= z_max

       G4double Pz = std::sqrt(Pz2);
       G4double zMin = (std::sqrt(HadronMassT2+Pz2) - Pz)/std::sqrt(StringMT2);
       G4double zMax = (std::sqrt(HadronMassT2+Pz2) + Pz)/std::sqrt(StringMT2);

       if (zMin >= zMax) return 0;		// have to start all over!
	
       G4double z = GetLightConeZ(zMin, zMax,
		                  string->GetDecayParton()->GetPDGEncoding(), pHadron,
		                  HadronPt.x(), HadronPt.y());      

       //... now compute hadron longitudinal momentum and energy
       // longitudinal hadron momentum component HadronPz

       HadronPt.setZ( 0.5* string->GetDecayDirection() *
		      (z * string->LightConeDecay() - 
		       HadronMassT2/(z * string->LightConeDecay())) );
       G4double HadronE  = 0.5* (z * string->LightConeDecay() + 
				 HadronMassT2/(z * string->LightConeDecay()) );

       G4LorentzVector * a4Momentum= new G4LorentzVector(HadronPt,HadronE);

       #ifdef debug_QGSMfragmentation
       G4cout<<"string->GetDecayDirection() string->LightConeDecay() "
             <<string->GetDecayDirection()<<" "<<string->LightConeDecay()<<G4endl;
       G4cout<<"HadronPt,HadronE "<<HadronPt<<" "<<HadronE<<G4endl;
       G4cout<<"Out of QGSM SplitEandP "<<G4endl;
       #endif

       return a4Momentum;
}

//----------------------------------------------------------------------------------------------------------

G4double G4QGSMFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int /* PartonEncoding */ ,  
                                            G4ParticleDefinition* /* pHadron */, G4double ptx , G4double pty)
{    
  G4double lambda = 2.0*(sqr(ptx)+sqr(pty))/sqr(GeV);

  #ifdef debug_QGSMfragmentation
  G4cout<<"GetLightConeZ zmin zmax Parton pHadron "<<zmin<<" "<<zmax<<" "/*<< PartonEncoding */
        <<" "/*<< pHadron->GetParticleName() */ <<G4endl;
  #endif

  G4double z(0.);    
  G4int DiQold(0), DiQnew(0);
  G4double d1(-1.0), d2(0.);
  G4double invD1(0.),invD2(0.), r1(0.),r2(0.),r12(0.);

  G4int absDecayQuarkCode = std::abs( DecayQuark );
  G4int absNewQuarkCode   = std::abs( NewQuark   );

  G4int q1(0), q2(0);
  //  q1 = absDecayQuarkCode/1000; q2 = (absDecayQuarkCode % 1000)/100;

  G4int qA(0), qB(0);
  //  qA = absNewQuarkCode/1000;   qB = (absNewQuarkCode % 1000)/100;

  if ( (absDecayQuarkCode < 6) && (absNewQuarkCode < 6) ) {
    d1 = FFq2q[absDecayQuarkCode-1][absNewQuarkCode-1][0]; d2 = FFq2q[absDecayQuarkCode-1][absNewQuarkCode-1][1];
  }

  if ( (absDecayQuarkCode < 6) && (absNewQuarkCode > 6) ) {
   qA = absNewQuarkCode/1000;   qB = (absNewQuarkCode % 1000)/100;   DiQnew = IndexDiQ[qA-1][qB-1];
   d1 = FFq2qq[absDecayQuarkCode-1][DiQnew][0]; d2 = FFq2q[absDecayQuarkCode-1][DiQnew][1];
  }

  if ( (absDecayQuarkCode > 6) && (absNewQuarkCode < 6) ) {
    q1 = absDecayQuarkCode/1000; q2 = (absDecayQuarkCode % 1000)/100; DiQold = IndexDiQ[q1-1][q2-1];
    d1 = FFqq2q[DiQold][absNewQuarkCode-1][0]; d2 = FFqq2q[DiQold][absNewQuarkCode-1][1];
  }

  if ( d1 < 0. ) {
    q1 = absDecayQuarkCode/1000; q2 = (absDecayQuarkCode % 1000)/100; DiQold = IndexDiQ[q1-1][q2-1];
    qA = absNewQuarkCode/1000;   qB = (absNewQuarkCode % 1000)/100;   DiQnew = IndexDiQ[qA-1][qB-1];
    d1 = FFqq2qq[DiQold][DiQnew][0]; d2 = FFqq2qq[DiQold][DiQnew][1];
  }

  d2 +=lambda;
  d1+=1.0; d2+=1.0;

  invD1=1./d1; invD2=1./d2;

  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = 0;
  do  // Jong's algorithm
  {
    r1=G4Pow::GetInstance()->powA(G4UniformRand(),invD1);
    r2=G4Pow::GetInstance()->powA(G4UniformRand(),invD2);
    r12=r1+r2;
    z=r1/r12;
  } while ( ( (r12 > 1.0) || !((zmin <= z)&&(z <= zmax))) && 
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */

  if ( loopCounter >= maxNumberOfLoops ) {
    z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
  }

  return z;
}

//-----------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::SplitLast(G4FragmentingString * string,
			              G4KineticTrackVector * LeftVector,
    				      G4KineticTrackVector * RightVector)
{
    //... perform last cluster decay

    G4ThreeVector ClusterVel =string->Get4Momentum().boostVector();
    G4double ResidualMass    =string->Mass();

    #ifdef debug_QGSMfragmentation
    G4cout<<"Split last-----------------------------------------"<<G4endl;
    G4cout<<"StrMass "<<ResidualMass<<" q's "
          <<string->GetLeftParton()->GetParticleName()<<" "
          <<string->GetRightParton()->GetParticleName()<<G4endl;
    #endif

    G4int cClusterInterrupt = 0;
    G4ParticleDefinition *LeftHadron = nullptr;
    G4ParticleDefinition *RightHadron = nullptr;
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;

    G4double LeftHadronMass(0.); G4double RightHadronMass(0.);
    do
    {
        if (cClusterInterrupt++ >= ClusterLoopInterrupt) return false;
        LeftHadronMass = -MaxMass; RightHadronMass = -MaxMass;

	G4ParticleDefinition * quark = nullptr;
	string->SetLeftPartonStable(); // to query quark contents..

	if (string->DecayIsQuark() && string->StableIsQuark() ) 
	{
	  //... there are quarks on cluster ends

	  G4int IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? -1 : +1; 
                // if we have a quark, we need antiquark or diquark

	  pDefPair QuarkPair = CreatePartonPair(IsParticle);
	  quark = QuarkPair.second;

	  LeftHadron= hadronizer->Build(QuarkPair.first, string->GetLeftParton());
          if ( LeftHadron == NULL ) continue;
          RightHadron = hadronizer->Build(string->GetRightParton(), quark);
          if ( RightHadron == NULL ) continue;    
	} else if( (!string->DecayIsQuark() &&  string->StableIsQuark() ) ||   
	           ( string->DecayIsQuark() && !string->StableIsQuark() )   ) {
	  //... there is a Diquark on one of cluster ends
	  G4int IsParticle;
	  if ( string->StableIsQuark() ) {
	    IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? -1 : +1; 
	  } else {
	    IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? +1 : -1;
	  }

          //G4double ProbSaS   = 1.0 - 2.0 * GetStrangeSuppress();
          //G4double ActualProb = ProbSaS * 1.4;
          //SetStrangenessSuppression((1.0-ActualProb)/2.0);

      	  pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
          //SetStrangenessSuppression((1.0-ProbSaS)/2.0);
      	  quark = QuarkPair.second;
      	  LeftHadron=hadronizer->Build(QuarkPair.first, string->GetLeftParton());
          if ( LeftHadron == NULL ) continue;
          RightHadron = hadronizer->Build(string->GetRightParton(), quark);
          if ( RightHadron == NULL ) continue;
        } else {  // Diquark and anti-diquark are on the string ends 
          if (cClusterInterrupt++ >= ClusterLoopInterrupt) return false;
          G4int LeftQuark1= string->GetLeftParton()->GetPDGEncoding()/1000;
          G4int LeftQuark2=(string->GetLeftParton()->GetPDGEncoding()/100)%10;
          G4int RightQuark1= string->GetRightParton()->GetPDGEncoding()/1000;
          G4int RightQuark2=(string->GetRightParton()->GetPDGEncoding()/100)%10;
          if (G4UniformRand()<0.5) {
            LeftHadron  =hadronizer->Build(FindParticle( LeftQuark1), FindParticle(RightQuark1));
            RightHadron =hadronizer->Build(FindParticle( LeftQuark2), FindParticle(RightQuark2));
          } else {
            LeftHadron  =hadronizer->Build(FindParticle( LeftQuark1), FindParticle(RightQuark2));
            RightHadron =hadronizer->Build(FindParticle( LeftQuark2), FindParticle(RightQuark1));
          }
	  if ( (LeftHadron == NULL) || (RightHadron == NULL) ) continue;
        }
        LeftHadronMass  = LeftHadron->GetPDGMass();
        RightHadronMass = RightHadron->GetPDGMass();
        //... repeat procedure, if mass of cluster is too low to produce hadrons
    } while ( ( ResidualMass <= LeftHadronMass + RightHadronMass )
              && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */

    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }

    //... compute hadron momenta and energies   
    G4LorentzVector  LeftMom, RightMom;
    G4ThreeVector    Pos;
    Sample4Momentum(&LeftMom , LeftHadron->GetPDGMass() , 
                    &RightMom, RightHadron->GetPDGMass(), ResidualMass);
    LeftMom.boost(ClusterVel);
    RightMom.boost(ClusterVel);

    #ifdef debug_QGSMfragmentation
    G4cout<<LeftHadron->GetParticleName()<<" "<<RightHadron->GetParticleName()<<G4endl;
    G4cout<<"Left  Hadrom P M "<<LeftMom<<" "<<LeftMom.mag()<<G4endl;
    G4cout<<"Right Hadrom P M "<<RightMom<<" "<<RightMom.mag()<<G4endl;
    #endif

    LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
    RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

    return true;
}

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::Sample4Momentum(G4LorentzVector* Mom    , G4double Mass    , 
                                          G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass) 
{
    #ifdef debug_QGSMfragmentation
    G4cout<<"Sample4Momentum Last-----------------------------------------"<<G4endl;
    G4cout<<"  StrMass "<<InitialMass<<" Mass1 "<<Mass<<" Mass2 "<<AntiMass<<G4endl;
    G4cout<<"  SumMass "<<Mass+AntiMass<<G4endl;
    #endif

    G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) - sqr(2.*Mass*AntiMass);
    G4double Pabs = (r_val > 0.)? std::sqrt(r_val)/(2.*InitialMass) : 0;

    //... sample unit vector       
    G4double pz = 1. - 2.*G4UniformRand();  
    G4double st     = std::sqrt(1. - pz * pz)*Pabs;
    G4double phi    = 2.*pi*G4UniformRand();
    G4double px = st*std::cos(phi);
    G4double py = st*std::sin(phi);
    pz *= Pabs;
    
    Mom->setPx(px); Mom->setPy(py); Mom->setPz(pz);
    Mom->setE(std::sqrt(Pabs*Pabs + Mass*Mass));

    AntiMom->setPx(-px); AntiMom->setPy(-py); AntiMom->setPz(-pz);
    AntiMom->setE (std::sqrt(Pabs*Pabs + AntiMass*AntiMass));
}

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::SetFFq2q()  // q-> q' + Meson (q anti q')
{
  for (G4int i=0; i < 5; i++) {
    FFq2q[i][0][0] = 2.0 ; FFq2q[i][0][1] = -arho + alft;  // q->d + (q dbar) Pi0, Eta, Eta', Rho0, omega
    FFq2q[i][1][0] = 2.0 ; FFq2q[i][1][1] = -arho + alft;  // q->u + (q ubar) Pi-, Rho-
    FFq2q[i][2][0] = 1.0 ; FFq2q[i][2][1] = -aphi + alft;  // q->s + (q sbar) K0, K*0
    FFq2q[i][3][0] = 1.0 ; FFq2q[i][3][1] = -aJPs + alft;  // q->c + (q+cbar) D-, D*-
    FFq2q[i][4][0] = 1.0 ; FFq2q[i][4][1] = -aUps + alft;  // q->b + (q bbar) EtaB, Upsilon
  }
}

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::SetFFq2qq()  // q-> anti (q1'q2') + Baryon (q + q1 + q2)
{
  for (G4int i=0; i < 5; i++) {
    FFq2qq[i][ 0][0] = 0.0 ; FFq2qq[i][ 0][1] = arho - 2.0*an    + alft;  // q->dd bar + (q dd)
    FFq2qq[i][ 1][0] = 0.0 ; FFq2qq[i][ 1][1] = arho - 2.0*an    + alft;  // q->ud bar + (q ud)
    FFq2qq[i][ 2][0] = 0.0 ; FFq2qq[i][ 2][1] = arho - 2.0*ala   + alft;  // q->sd bar + (q sd) 
    FFq2qq[i][ 3][0] = 0.0 ; FFq2qq[i][ 3][1] = arho - 2.0*alaC  + alft;  // q->cd bar + (q cd)
    FFq2qq[i][ 4][0] = 0.0 ; FFq2qq[i][ 4][1] = arho - 2.0*alaB  + alft;  // q->bd bar + (q bd)
    FFq2qq[i][ 5][0] = 0.0 ; FFq2qq[i][ 5][1] = arho - 2.0*an    + alft;  // q->uu bar + (q uu)
    FFq2qq[i][ 6][0] = 0.0 ; FFq2qq[i][ 6][1] = arho - 2.0*ala   + alft;  // q->su bar + (q su)
    FFq2qq[i][ 7][0] = 0.0 ; FFq2qq[i][ 7][1] = arho - 2.0*alaC  + alft;  // q->cu bar + (q cu)
    FFq2qq[i][ 8][0] = 0.0 ; FFq2qq[i][ 8][1] = arho - 2.0*alaB  + alft;  // q->bu bar + (q bu)
    FFq2qq[i][ 9][0] = 0.0 ; FFq2qq[i][ 9][1] = arho - 2.0*aXi   + alft;  // q->ss bar + (q ss)    
    FFq2qq[i][10][0] = 0.0 ; FFq2qq[i][10][1] = arho - 2.0*aXiC  + alft;  // q->cs bar + (q cs)
    FFq2qq[i][11][0] = 0.0 ; FFq2qq[i][11][1] = arho - 2.0*aXiB  + alft;  // q->bs bar + (q bc)
    FFq2qq[i][12][0] = 0.0 ; FFq2qq[i][12][1] = arho - 2.0*aXiCC + alft;  // q->cc bar + (q cc) 
    FFq2qq[i][13][0] = 0.0 ; FFq2qq[i][13][1] = arho - 2.0*aXiCB + alft;  // q->bc bar + (q bc)
    FFq2qq[i][14][0] = 0.0 ; FFq2qq[i][14][1] = arho - 2.0*aXiBB + alft;  // q->bb bar + (q bb)
  }
}

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::SetFFqq2q()  // q1q2-> anti(q') + Baryon (q1 + q2 + q')
{
  for (G4int i=0; i < 15; i++) {   
    FFqq2q[i][0][0] = 2.0*(arho - an); FFqq2q[i][0][1] = -arho + alft;
    FFqq2q[i][1][0] = 2.0*(arho - an); FFqq2q[i][1][1] = -arho + alft;
    FFqq2q[i][2][0] = 2.0*(arho - an); FFqq2q[i][2][1] = -aphi + alft;
    FFqq2q[i][3][0] = 2.0*(arho - an); FFqq2q[i][3][1] = -aJPs + alft;
    FFqq2q[i][4][0] = 2.0*(arho - an); FFqq2q[i][4][1] = -aUps + alft;
  }
}

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::SetFFqq2qq()  // q1(q2)-> q'(q2) + Meson(q1 anti q')
{
  for (G4int i=0; i < 15; i++) {
    FFqq2qq[i][0][0] = 0.  ;  FFqq2qq[i][0][1] = 2.0*arho - 2.0*an -arho + alft;  // dd -> dd + Pi0 (d d bar)
    FFqq2qq[i][1][0] = 0.  ;  FFqq2qq[i][1][1] = 2.0*arho - 2.0*an -arho + alft;  // dd -> ud + Pi- (d u bar)
    FFqq2qq[i][2][0] = 0.  ;  FFqq2qq[i][2][1] = 2.0*arho - 2.0*an -aphi + alft;  // dd -> sd + K0  (d s bar)
    FFqq2qq[i][3][0] = 0.  ;  FFqq2qq[i][3][1] = 2.0*arho - 2.0*an -aJPs + alft;  // dd -> cd + D-  (d c bar)
    FFqq2qq[i][4][0] = 0.  ;  FFqq2qq[i][4][1] = 2.0*arho - 2.0*an -aUps + alft;  // dd -> bd + B0  (d b bar)
  }
}

