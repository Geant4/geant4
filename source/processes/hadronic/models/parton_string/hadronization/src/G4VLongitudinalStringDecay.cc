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
// $Id: G4VLongitudinalStringDecay.cc,v 1.5 2006/06/29 20:55:09 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
//               redesign  Gunter Folger, August/September 2001
// -----------------------------------------------------------------------------
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4VLongitudinalStringDecay.hh"
#include "G4FragmentingString.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleChange.hh"
#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4Gluons.hh"

//********************************************************************************
// Constructors

G4VLongitudinalStringDecay::G4VLongitudinalStringDecay()
{
   MassCut  = 0.35*GeV; 
   ClusterMass = 0.15*GeV;

   SmoothParam      = 0.9; 
   StringLoopInterrupt    = 1000;
   ClusterLoopInterrupt   =  500;

// Changable Parameters below.
   
   SigmaQT = 0.5 * GeV;
   
   StrangeSuppress  = 0.44;    //  27 % strange quarks produced, ie. u:d:s=1:1:0.27
   DiquarkSuppress  = 0.1;
   DiquarkBreakProb = 0.1;
   
   //... pspin_meson is probability to create vector meson 
   pspin_meson = 0.5;

   //... pspin_barion is probability to create 3/2 barion 
   pspin_barion = 0.5;

   //... vectorMesonMix[] is quark mixing parameters for vector mesons (Variable spin = 3)
   vectorMesonMix.resize(6);
   vectorMesonMix[0] = 0.5;
   vectorMesonMix[1] = 0.0;
   vectorMesonMix[2] = 0.5;
   vectorMesonMix[3] = 0.0;
   vectorMesonMix[4] = 1.0;
   vectorMesonMix[5] = 1.0; 

   //... scalarMesonMix[] is quark mixing parameters for scalar mesons (Variable spin=1)
   scalarMesonMix.resize(6);
   scalarMesonMix[0] = 0.5; 
   scalarMesonMix[1] = 0.25; 
   scalarMesonMix[2] = 0.5; 
   scalarMesonMix[3] = 0.25; 
   scalarMesonMix[4] = 1.0; 
   scalarMesonMix[5] = 0.5; 

// Parameters may be changed until the first fragmentation starts
   PastInitPhase=false;
   hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
}
   

G4VLongitudinalStringDecay::~G4VLongitudinalStringDecay()
   {
   delete hadronizer;
   }

//=============================================================================================-------------

// Operators

//const  & G4VLongitudinalStringDecay::operator=(const G4VLongitudinalStringDecay &)
//    {
//    }

//----------------------------------------------------------------------------------------------------------

int G4VLongitudinalStringDecay::operator==(const G4VLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::operator== forbidden");
	return false;
    }

//----------------------------------------------------------------------------------------------------------

int G4VLongitudinalStringDecay::operator!=(const G4VLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::operator!= forbidden");
	return true;
    }

//==========================================================================================================

G4int G4VLongitudinalStringDecay::SampleQuarkFlavor(void)
   {
   return (1 + (int)(G4UniformRand()/StrangeSuppress));
   }

//----------------------------------------------------------------------------------------------------------

G4VLongitudinalStringDecay::pDefPair G4VLongitudinalStringDecay::CreatePartonPair(G4int NeedParticle,G4bool AllowDiquarks)
{
//  NeedParticle = +1 for Particle, -1 for Antiparticle

    if ( AllowDiquarks && G4UniformRand() < DiquarkSuppress )
    {
      // Create a Diquark - AntiDiquark pair , first in pair is anti to IsParticle
      G4int q1  = SampleQuarkFlavor();
      G4int q2  = SampleQuarkFlavor();
      G4int spin = (q1 != q2 && G4UniformRand() <= 0.5)? 1 : 3;
                                     //   convention: quark with higher PDG number is first
      G4int PDGcode = (std::max(q1,q2) * 1000 + std::min(q1,q2) * 100 + spin) * NeedParticle;
      return pDefPair (FindParticle(-PDGcode),FindParticle(PDGcode));
      

    } else {
      // Create a Quark - AntiQuark pair, first in pair  IsParticle
      G4int PDGcode=SampleQuarkFlavor()*NeedParticle;
      return pDefPair (FindParticle(PDGcode),FindParticle(-PDGcode));
    }

}

//----------------------------------------------------------------------------------------------------------

G4ThreeVector G4VLongitudinalStringDecay::SampleQuarkPt()
   {
   G4double Pt = -std::log(G4UniformRand());
   Pt = SigmaQT * std::sqrt(Pt);
   G4double phi = 2.*pi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
   }

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector* Hadrons)
   {
   // `yo-yo` formation time
   const G4double kappa = 1.0 * GeV/fermi;
   for(size_t c1 = 0; c1 < Hadrons->size(); c1++)
      {
      G4double SumPz = 0; 
      G4double SumE  = 0;
      for(size_t c2 = 0; c2 < c1; c2++)
         {
         SumPz += Hadrons->operator[](c2)->Get4Momentum().pz();
         SumE  += Hadrons->operator[](c2)->Get4Momentum().e();   
         } 
      G4double HadronE  = Hadrons->operator[](c1)->Get4Momentum().e();
      G4double HadronPz = Hadrons->operator[](c1)->Get4Momentum().pz();
      Hadrons->operator[](c1)->SetFormationTime((theInitialStringMass - 2.*SumPz + HadronE - HadronPz)/(2.*kappa));
      G4ThreeVector aPosition(0, 0,     (theInitialStringMass - 2.*SumE  - HadronE + HadronPz)/(2.*kappa));
      Hadrons->operator[](c1)->SetPosition(aPosition);
      } 
   }

//----------------------------------------------------------------------------------------------------------

/*
void G4VLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector* Hadrons)
   {
   // 'constituent' formation time 
   const G4double kappa = 1.0 * GeV/fermi;
   for(G4int c1 = 0; c1 < Hadrons->length(); c1++)
      {
      G4double SumPz = 0; 
      G4double SumE  = 0;
      for(G4int c2 = 0; c2 <= c1; c2++)
         {
         SumPz += Hadrons->at(c2)->Get4Momentum().pz();
         SumE  += Hadrons->at(c2)->Get4Momentum().e();   
         } 
      Hadrons->at(c1)->SetFormationTime((theInitialStringMass - 2.*SumPz)/(2.*kappa));
      G4ThreeVector aPosition(0, 0,     (theInitialStringMass - 2.*SumE)/(2.*kappa));
      Hadrons->at(c1)->SetPosition(aPosition);
      } 
   c1 = Hadrons->length()-1;   
   Hadrons->at(c1)->SetFormationTime(Hadrons->at(c1-1)->GetFormationTime());
   Hadrons->at(c1)->SetPosition(Hadrons->at(c1-1)->GetPosition());
   }
*/

//----------------------------------------------------------------------------------------------------------

G4ParticleDefinition *
		G4VLongitudinalStringDecay::QuarkSplitup(G4ParticleDefinition*
		decay, G4ParticleDefinition *&created)
{
    G4int IsParticle=(decay->GetPDGEncoding()>0) ? -1 : +1; // if we have a quark, we need antiquark (or diquark)
    pDefPair QuarkPair = CreatePartonPair(IsParticle);
    created = QuarkPair.second;
    return hadronizer->Build(QuarkPair.first, decay);
    
}

//----------------------------------------------------------------------------------------------------------

G4ParticleDefinition *G4VLongitudinalStringDecay::DiQuarkSplitup(
							G4ParticleDefinition* decay,
							G4ParticleDefinition *&created)
{
   //... can Diquark break or not? 
   if (G4UniformRand() < DiquarkBreakProb ){
   //... Diquark break

      G4int stableQuarkEncoding = decay->GetPDGEncoding()/1000;
      G4int decayQuarkEncoding = (decay->GetPDGEncoding()/100)%10;
      if (G4UniformRand() < 0.5)
         {
         G4int Swap = stableQuarkEncoding;
         stableQuarkEncoding = decayQuarkEncoding;
         decayQuarkEncoding = Swap;
         }

      G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1; 
			// if we have a quark, we need antiquark)
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      //... Build new Diquark
      G4int QuarkEncoding=QuarkPair.second->GetPDGEncoding();
      G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
      G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
      created = FindParticle(NewDecayEncoding);
      G4ParticleDefinition * decayQuark=FindParticle(decayQuarkEncoding);
      
      return hadronizer->Build(QuarkPair.first, decayQuark);
   
   } else {
   //... Diquark does not break
 
      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1; 
			// if we have a diquark, we need quark)
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      created = QuarkPair.second;

      return hadronizer->Build(QuarkPair.first, decay);
   }
}

//-----------------------------------------------------------------------------------------

G4KineticTrack * G4VLongitudinalStringDecay::Splitup(
		        G4FragmentingString *string, 
			G4FragmentingString *&newString)
{

       //... random choice of string end to use for creating the hadron (decay)   
       SideOfDecay = (G4UniformRand() < 0.5)? 1: -1;
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
       	   HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);
       } else {
           HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
       }      
// create new String from old, ie. keep Left and Right order, but replace decay
       G4LorentzVector* HadronMomentum=SplitEandP(HadronDefinition, string);
	
       G4KineticTrack * Hadron =0;
       if ( HadronMomentum != 0 ) {    

	   G4ThreeVector   Pos;
	   Hadron = new G4KineticTrack(HadronDefinition, 0,Pos, *HadronMomentum);
 
	   newString=new G4FragmentingString(*string,newStringEnd,
	   				HadronMomentum);
	   
	   delete HadronMomentum;
       }      
       return Hadron;
}

//-----------------------------------------------------------------------------------------

G4LorentzVector * G4VLongitudinalStringDecay::SplitEandP(G4ParticleDefinition * pHadron,
	G4FragmentingString * string)
{
       G4double HadronMass = pHadron->GetPDGMass();

       // calculate and assign hadron transverse momentum component HadronPx andHadronPy
       G4ThreeVector thePt;
       thePt=SampleQuarkPt();
       G4ThreeVector HadronPt = thePt +string->DecayPt();
       HadronPt.setZ(0);
       //...  sample z to define hadron longitudinal momentum and energy
       //... but first check the available phase space
       G4double DecayQuarkMass2  = sqr(string->GetDecayParton()->GetPDGMass());
       G4double HadronMass2T = sqr(HadronMass) + HadronPt.mag2();
       if (DecayQuarkMass2 + HadronMass2T >= SmoothParam*(string->Mass2()) ) 
          return 0;		// have to start all over!

       //... then compute allowed z region  z_min <= z <= z_max 
 
       G4double zMin = HadronMass2T/(string->Mass2());
       G4double zMax = 1. - DecayQuarkMass2/(string->Mass2());
       if (zMin >= zMax) return 0;		// have to start all over!
	
       G4double z = GetLightConeZ(zMin, zMax,
		       string->GetDecayParton()->GetPDGEncoding(), pHadron,
		       HadronPt.x(), HadronPt.y());      
       
       //... now compute hadron longitudinal momentum and energy
       // longitudinal hadron momentum component HadronPz

        HadronPt.setZ(0.5* string->GetDecayDirection() *
			(z * string->LightConeDecay() - 
			 HadronMass2T/(z * string->LightConeDecay())));
        G4double HadronE  = 0.5* (z * string->LightConeDecay() + 
				  HadronMass2T/(z * string->LightConeDecay()));

       G4LorentzVector * a4Momentum= new G4LorentzVector(HadronPt,HadronE);

       return a4Momentum;
}


//-----------------------------------------------------------------------------------------

G4bool G4VLongitudinalStringDecay::SplitLast(G4FragmentingString * string,
					     G4KineticTrackVector * LeftVector,
    					     G4KineticTrackVector * RightVector)
{
    //... perform last cluster decay
    G4ThreeVector ClusterVel =string->Get4Momentum().boostVector();
    G4double ResidualMass    =string->Mass(); 
    G4double ClusterMassCut = ClusterMass;
    G4int cClusterInterrupt = 0;
    G4ParticleDefinition * LeftHadron, * RightHadron;
    do
    {
        if (cClusterInterrupt++ >= ClusterLoopInterrupt)
        {
          return false;
        }
	G4ParticleDefinition * quark = NULL;
	string->SetLeftPartonStable(); // to query quark contents..
	if (string->DecayIsQuark() && string->StableIsQuark() ) 
	{
	   //... there are quarks on cluster ends
		LeftHadron= QuarkSplitup(string->GetLeftParton(), quark);
	} else {
	   //... there is a Diquark on cluster ends
		G4int IsParticle;
		if ( string->StableIsQuark() ) {
		  IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? -1 : +1; 
		} else {
		  IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? +1 : -1;
		}
      		pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      		quark = QuarkPair.second;
      		LeftHadron=hadronizer->Build(QuarkPair.first, string->GetLeftParton());
	}
        RightHadron = hadronizer->Build(string->GetRightParton(), quark);

       //... repeat procedure, if mass of cluster is too low to produce hadrons
       //... ClusterMassCut = 0.15*GeV model parameter
	if ( quark->GetParticleSubType()== "quark" ) {ClusterMassCut = 0.;}
	else {ClusterMassCut = ClusterMass;}
    } 
    while (ResidualMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()  + ClusterMassCut);

    //... compute hadron momenta and energies   
    G4LorentzVector  LeftMom, RightMom;
    G4ThreeVector    Pos;
    Sample4Momentum(&LeftMom, LeftHadron->GetPDGMass(), &RightMom, RightHadron->GetPDGMass(), ResidualMass);
    LeftMom.boost(ClusterVel);
    RightMom.boost(ClusterVel);
    LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
    RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

    return true;

}

//----------------------------------------------------------------------------------------------------------

G4KineticTrackVector* G4VLongitudinalStringDecay::FragmentString(const G4ExcitedString& theString)
{
//    Can no longer modify Parameters for Fragmentation.
	PastInitPhase=true;
	
// 	check if string has enough mass to fragment...
	G4KineticTrackVector * LeftVector=LightFragmentationTest(&theString);
	if ( LeftVector != 0 ) return LeftVector;
	
	LeftVector = new G4KineticTrackVector;
	G4KineticTrackVector * RightVector=new G4KineticTrackVector;

// this should work but its only a semi deep copy. %GF	G4ExcitedString theStringInCMS(theString);
        G4ExcitedString *theStringInCMS=CPExcited(theString);
	G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();

	G4bool success=false, inner_sucess=true;
	G4int attempt=0;
	while ( !success && attempt++ < StringLoopInterrupt )
	{
		G4FragmentingString *currentString=new G4FragmentingString(*theStringInCMS);

		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		RightVector->clear();
		
		inner_sucess=true;  // set false on failure..
		while (! StopFragmenting(currentString) )
		{  // Split current string into hadron + new string
			G4FragmentingString *newString=0;  // used as output from SplitUp...
			G4KineticTrack * Hadron=Splitup(currentString,newString);
			if ( Hadron != 0 && IsFragmentable(newString)) 
			{
			   if ( currentString->GetDecayDirection() > 0 )
				   LeftVector->push_back(Hadron);
       			   else
	  			   RightVector->push_back(Hadron);
			   delete currentString;
			   currentString=newString;
			} else {
			 // abandon ... start from the beginning
			   if (newString) delete newString;
			   if (Hadron)    delete Hadron;
			   inner_sucess=false;
			   break;
			}
		} 
		// Split current string into 2 final Hadrons
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
	while(!RightVector->empty())
	{
	    LeftVector->push_back(RightVector->back());
	    RightVector->erase(RightVector->end()-1);
	}
	delete RightVector;

	CalculateHadronTimePosition(theString.Get4Momentum().mag(), LeftVector);

	G4LorentzRotation toObserverFrame(toCms.inverse());

	for(size_t C1 = 0; C1 < LeftVector->size(); C1++)
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

G4ExcitedString *G4VLongitudinalStringDecay::CPExcited(const G4ExcitedString & in)
{
	G4Parton *Left=new G4Parton(*in.GetLeftParton());
	G4Parton *Right=new G4Parton(*in.GetRightParton());
	return new G4ExcitedString(Left,Right,in.GetDirection());
}

G4double G4VLongitudinalStringDecay::FragmentationMass(
		const G4FragmentingString *
		const string,
		Pcreate build,
		pDefPair * pdefs)
{
	
        G4double mass;

	if ( build==0 ) build=&G4HadronBuilder::BuildLowSpin;

        G4ParticleDefinition *Hadron1, *Hadron2=0;

        if (!string->FourQuarkString() )
        {
           // spin 0 meson or spin 1/2 barion will be built

           Hadron1 = (hadronizer->*build)(string->GetLeftParton(),
			              string->GetRightParton());
           mass= (Hadron1)->GetPDGMass();
        } else
        {
           //... string is qq--qqbar: Build two stable hadrons,
           //...    with extra uubar or ddbar quark pair
	   G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
	   if (string->GetLeftParton()->GetPDGEncoding() < 0) iflc = -iflc;

	   //... theSpin = 4; spin 3/2 baryons will be built
	   Hadron1 = (hadronizer->*build)(string->GetLeftParton(),FindParticle(iflc));
	   Hadron2 =(hadronizer->*build)(string->GetRightParton(),FindParticle(-iflc));
           mass = (Hadron1)->GetPDGMass() + (Hadron2)->GetPDGMass();
        }
	
	if ( pdefs != 0 ) 
	{ // need to return hadrons as well....
	   pdefs->first  = Hadron1;
	   pdefs->second = Hadron2;
	}
	   
        return mass;
}

//----------------------------------------------------------------------------------------------------------

G4bool G4VLongitudinalStringDecay::IsFragmentable(const G4FragmentingString * const string)
{
	return sqr(FragmentationMass(string)+MassCut) <
			string->Mass2();
}

//----------------------------------------------------------------------------------------------------------

G4bool G4VLongitudinalStringDecay::StopFragmenting(const G4FragmentingString * const string)
{
	return
         sqr(FragmentationMass(string,&G4HadronBuilder::BuildHighSpin)+MassCut) >
         string->Get4Momentum().mag2();
}

//----------------------------------------------------------------------------------------------------------

G4KineticTrackVector* G4VLongitudinalStringDecay::LightFragmentationTest(const
		G4ExcitedString * const string)
{
   // Check string decay threshold
		
	G4KineticTrackVector * result=0;  // return 0 when string exceeds the mass cut
	
	pDefPair hadrons((G4ParticleDefinition *)0,(G4ParticleDefinition *)0);
	G4FragmentingString aString(*string);
	if ( sqr(FragmentationMass(&aString,0,&hadrons)+MassCut) < aString.Mass2()) {
		return 0;
	}
	
	result=new G4KineticTrackVector;
        
	if ( hadrons.second ==0 )
	{
	      	 // Substitute string by light hadron, Note that Energy is not conserved here!

	       G4ThreeVector Mom3 = string->Get4Momentum().vect();
	       G4LorentzVector Mom(Mom3, 
	       			   std::sqrt(Mom3.mag2() + sqr(hadrons.first->GetPDGMass())));
               result->push_back(new G4KineticTrack(hadrons.first, 0, string->GetPosition(), Mom));
	} else 
	{
	   //... string was qq--qqbar type: Build two stable hadrons,
	       G4LorentzVector  Mom1, Mom2;
	       Sample4Momentum(&Mom1, hadrons.first->GetPDGMass(), 
			       &Mom2,hadrons.second->GetPDGMass(),
			        string->Get4Momentum().mag());
	       result->push_back(new G4KineticTrack(hadrons.first, 0, string->GetPosition(), Mom1));
	       result->push_back(new G4KineticTrack(hadrons.second, 0, string->GetPosition(), Mom2));
               G4ThreeVector Velocity = string->Get4Momentum().boostVector();
               result->Boost(Velocity);          
	}

	return result;
	
}

//----------------------------------------------------------------------------------------------------------

G4ParticleDefinition* G4VLongitudinalStringDecay::FindParticle(G4int Encoding) 
   {
   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);
      if (ptr == NULL)
       {
       G4cout << "Particle with encoding "<<Encoding<<" does not exist!!!"<<G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Check your particle table");
       }
   return ptr;    
   }

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass) 
    {
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

void G4VLongitudinalStringDecay::SetSigmaTransverseMomentum(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetSigmaTransverseMomentum after FragmentString() not allowed");
	} else {
		SigmaQT = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetStrangenessSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetStrangenessSuppression after FragmentString() not allowed");
	} else {
		StrangeSuppress = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetDiquarkSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetDiquarkSuppression after FragmentString() not allowed");
	} else {
		DiquarkSuppress = aValue;
	}
}

void G4VLongitudinalStringDecay::SetDiquarkBreakProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetDiquarkBreakProbability after FragmentString() not allowed");
	} else {
		DiquarkBreakProb = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetVectorMesonProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetVectorMesonProbability after FragmentString() not allowed");
	} else {
		pspin_meson = aValue;
		delete hadronizer;
		hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability after FragmentString() not allowed");
	} else {
		pspin_barion = aValue;
		delete hadronizer;
		hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetScalarMesonMixings(std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetScalarMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetScalarMesonMixings( argument Vector too small");
	  scalarMesonMix[0] = aVector[0];
	  scalarMesonMix[1] = aVector[1];
	  scalarMesonMix[2] = aVector[2];
	  scalarMesonMix[3] = aVector[3];
	  scalarMesonMix[4] = aVector[4];
	  scalarMesonMix[5] = aVector[5];
	  delete hadronizer;
	  hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetVectorMesonMixings(std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetVectorMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetVectorMesonMixings( argument Vector too small");
	  vectorMesonMix[0] = aVector[0];
	  vectorMesonMix[1] = aVector[1];
	  vectorMesonMix[2] = aVector[2];
	  vectorMesonMix[3] = aVector[3];
	  vectorMesonMix[4] = aVector[4];
	  vectorMesonMix[5] = aVector[5];
	  delete hadronizer;
	  hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
  
	}
}	
//*******************************************************************************************************
