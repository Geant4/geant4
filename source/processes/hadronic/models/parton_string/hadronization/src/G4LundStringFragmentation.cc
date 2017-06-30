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
// $Id: G4LundStringFragmentation.cc 103022 2017-03-09 11:12:29Z gcosmo $
// GEANT4 tag $Name:  $ 1.8
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4LundStringFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4FragmentingString.hh"
#include "G4DiQuarks.hh"
#include "G4Quarks.hh"

#include "G4Exp.hh"
#include "G4Pow.hh"

//#define debug_LUNDfragmentation

// Class G4LundStringFragmentation 
//*************************************************************************************

G4LundStringFragmentation::G4LundStringFragmentation()
{
    SetMassCut(160.*MeV);   // For LightFragmentationTest it is required
                            // that no one pi-meson can be produced.

    SigmaQT = 0.435 * GeV;  // Uzhi 28 May 2016 For mesons

    SetStringTensionParameter(1.*GeV/fermi);
    SetDiquarkBreakProbability(0.05); 
    SetStrangenessSuppression(0.450); //(0.473); 			// Uzhi June 2016 0.45 -> 0.473
    SetDiquarkSuppression(0.056);     			// Uzhi June 2016 0.06 -> 0.05

// For treating of small string decays
   SetMinMasses();
}

//--------------------------------------------------------------------------------------
G4KineticTrackVector* G4LundStringFragmentation::FragmentString(const G4ExcitedString& theString)
{
	// Can no longer modify Parameters for Fragmentation.
	PastInitPhase=true;

	G4FragmentingString  aString(theString);
	SetMinimalStringMass(&aString);

#ifdef debug_LUNDfragmentation
  G4cout<<G4endl<<"LUND StringFragmentation ------------------------------------"<<G4endl;
  G4cout<<G4endl<<"LUND StringFragm: String Mass "
                             <<theString.Get4Momentum().mag()<<G4endl
                             <<" 4Mom "<<theString.Get4Momentum()
                             <<"------------------------------------"<<G4endl;
  G4cout<<"String ends Direct "<<theString.GetLeftParton()->GetPDGcode()<<" "
                               <<theString.GetRightParton()->GetPDGcode()<<" "
                               <<theString.GetDirection()<< G4endl;
  G4cout<<"Left  mom "<<theString.GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"Right mom "<<theString.GetRightParton()->Get4Momentum()<<G4endl;
  G4cout<<"Check for Fragmentation "<<G4endl;
#endif

	G4KineticTrackVector * LeftVector(0);

	if(!aString.FourQuarkString() && !IsFragmentable(&aString))
	{
#ifdef debug_LUNDfragmentation
  G4cout<<"Non fragmentable - the string is converted to one hadron "<<G4endl;
#endif
//	SetMassCut(160.*MeV); // For LightFragmentationTest it is required
//	                      // that no one pi-meson can be produced.

		G4double Mcut = GetMassCut();
		SetMassCut(10000.*MeV);
		LeftVector=LightFragmentationTest(&theString);
		SetMassCut(Mcut);

		LeftVector->operator[](0)->SetFormationTime(theString.GetTimeOfCreation());
		LeftVector->operator[](0)->SetPosition(theString.GetPosition());

		if(LeftVector->size() > 1)
                {
		        // 2 hadrons created from qq-qqbar are stored
			LeftVector->operator[](1)->SetFormationTime(theString.GetTimeOfCreation());
			LeftVector->operator[](1)->SetPosition(theString.GetPosition());
		}
		return LeftVector;
	}  // end of if(!IsFragmentable(&aString))

#ifdef debug_LUNDfragmentation
  G4cout<<"The string will be fragmented. "<<G4endl;
#endif

	// The string can fragment. At least two particles can be produced.
			       LeftVector =new G4KineticTrackVector;
	G4KineticTrackVector * RightVector=new G4KineticTrackVector;

/* Uzhi Dec. 2016
	G4ExcitedString *theStringInCMS=CopyExcited(theString);

#ifdef debug_LUNDfragmentation
  G4cout<<"CMS Left  mom "<<theStringInCMS->GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"CMS Right mom "<<theStringInCMS->GetRightParton()->Get4Momentum()<<G4endl;
#endif

	G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();

#ifdef debug_LUNDfragmentation
  G4cout<<"aligCMS Left  mom "<<theStringInCMS->GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"aligCMS Right mom "<<theStringInCMS->GetRightParton()->Get4Momentum()<<G4endl;
  G4cout<<G4endl<<"LUND StringFragm: String Mass "
                             <<theStringInCMS->Get4Momentum().mag()<<" Pz Lab "
                             <<theStringInCMS->Get4Momentum().pz()
                             <<"------------------------------------"<<G4endl;
  G4cout<<"String ends and Direction "<<theStringInCMS->GetLeftParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetRightParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetDirection()<< G4endl;
#endif
*/ // Uzhi Dec. 2016

        G4bool success = Loop_toFragmentString(theString, LeftVector, RightVector);

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

	return LeftVector;
}

//----------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::IsFragmentable(const G4FragmentingString * const string)
{
	SetMinimalStringMass(string);
//  G4cout<<"MinM StrM "<<MinimalStringMass<<" "<< string->Get4Momentum().mag()<<G4endl;
	return MinimalStringMass < string->Get4Momentum().mag();
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::Loop_toFragmentString(
			const G4ExcitedString  &theString,
                        G4KineticTrackVector * & LeftVector,
                        G4KineticTrackVector * & RightVector )
{
#ifdef debug_LUNDfragmentation
  G4cout<<"Loop_toFrag "<<theString.GetLeftParton()->GetPDGcode()<<" "
                        <<theString.GetLeftParton()->Get4Momentum()<<G4endl
        <<"            "<<theString.GetRightParton()->GetPDGcode()<<" "
                        <<theString.GetRightParton()->Get4Momentum()<<G4endl
        <<"Direction   "<<theString.GetDirection()<< G4endl;
#endif

	G4bool final_success=false;
	G4bool inner_success=true;
	G4int attempt=0;
	while ( ! final_success && attempt++ < StringLoopInterrupt )
	{       // If the string fragmentation does not be happend, 
	        // repeat the fragmentation.

		G4FragmentingString *currentString=new G4FragmentingString(theString);
		G4LorentzRotation toCms, toObserverFrame;

		//G4cout<<"Main loop start whilecounter "<<attempt<<G4endl;

		// Cleaning up the previously produced hadrons
		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		RightVector->clear();

		// Main fragmentation loop until the string will not be able to fragment
		inner_success=true;  // set false on failure.
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = -1;

		while ( (! StopFragmenting(currentString)) && ++loopCounter < maxNumberOfLoops )
		{       // Split current string into hadron + new string
#ifdef debug_LUNDfragmentation
  G4cout<<"The string can fragment. "<<G4endl;;
//G4cout<<"1 "<<currentString->GetDecayDirection()<<G4endl;
#endif
			G4FragmentingString *newString=0;  // used as output from SplitUp.

			toCms=currentString->TransformToAlignedCms();
toObserverFrame= toCms.inverse(); // Uzhi Jan. 2017
#ifdef debug_LUNDfragmentation
  G4cout<<"CMS Left  mom "<<currentString->GetPleft()<<G4endl;
  G4cout<<"CMS Right mom "<<currentString->GetPright()<<G4endl;
  G4cout<<"CMS String M  "<<currentString->GetPstring()<<G4endl;
#endif

			G4KineticTrack * Hadron=Splitup(currentString,newString);
//G4cout<<"*Hadron --- "<<Hadron<<G4endl;
			if ( Hadron != 0 )  // Store the hadron                               
			{
#ifdef debug_LUNDfragmentation
  G4cout<<"Hadron prod at fragm. "<<Hadron->GetDefinition()->GetParticleName()<<G4endl;
//G4cout<<"2 "<<currentString->GetDecayDirection()<<G4endl;
#endif

//Uzhi Jan 2017				toObserverFrame= toCms.inverse();
				Hadron->Set4Momentum(toObserverFrame*Hadron->Get4Momentum());

				G4double TimeOftheStringCreation=theString.GetTimeOfCreation();
				G4ThreeVector PositionOftheStringCreation(theString.GetPosition());

				G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
				G4LorentzVector Momentum = toObserverFrame*Coordinate;
				Hadron->SetFormationTime(TimeOftheStringCreation + Momentum.e() - fermi/c_light);
				G4ThreeVector aPosition(Momentum.vect());
				Hadron->SetPosition(PositionOftheStringCreation+aPosition);


				if ( currentString->GetDecayDirection() > 0 )
                                {
					LeftVector->push_back(Hadron); 
                                } else
                                {
					RightVector->push_back(Hadron);
                                }
				delete currentString;
				currentString=newString;             // Uzhi June 2016
//Uzhi Jan. 2017				currentString->LorentzRotate(toObserverFrame);
			} else {
                          if ( newString ) delete newString;  //AR-09Mar2017 Fixed memory leak
                        }

currentString->LorentzRotate(toObserverFrame); // Uzhi Jan. 2017
		};   // while ( (! StopFragmenting(currentString)) ...
                if ( loopCounter >= maxNumberOfLoops ) {
                  inner_success=false;
                }
		// Split remaining string into 2 final hadrons.
#ifdef debug_LUNDfragmentation
  if(inner_success) G4cout<<"Split remaining string into 2 final hadrons."<<G4endl;
#endif

		if ( inner_success && SplitLast(currentString, LeftVector, RightVector) )
		{
		 final_success = true;
		}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++
	toObserverFrame=toCms.inverse();

	G4double TimeOftheStringCreation=theString.GetTimeOfCreation();
	G4ThreeVector PositionOftheStringCreation(theString.GetPosition());

	for(size_t C1 = 0; C1 < LeftVector->size(); C1++)
	{
		G4KineticTrack* Hadron = LeftVector->operator[](C1);
		G4LorentzVector Momentum = Hadron->Get4Momentum();
		//G4cout<<"Hadron "<<Hadron->GetDefinition()->GetParticleName()<<" "<<Momentum<<G4endl;
		Momentum = toObserverFrame*Momentum;
		Hadron->Set4Momentum(Momentum);

		G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
		Momentum = toObserverFrame*Coordinate;
		Hadron->SetFormationTime(TimeOftheStringCreation + Momentum.e() - fermi/c_light);
		G4ThreeVector aPosition(Momentum.vect());
		Hadron->SetPosition(PositionOftheStringCreation+aPosition);
	};
*/ //+++++++++++++++++++++++++++++++++=

		delete currentString;
	}  // End of the loop where we try to fragment the string.
	return final_success;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::StopFragmenting(const G4FragmentingString * const string)
{
	SetMinimalStringMass(string); 
	if (string->FourQuarkString())
	{
		return G4UniformRand() < G4Exp(-0.0005*(string->Mass() - MinimalStringMass));
	} else {
           
		G4bool Result = G4UniformRand() < 
				G4Exp(-0.66e-6*(string->Mass()*string->Mass() - MinimalStringMass*MinimalStringMass));
#ifdef debug_LUNDfragmentation 
  G4cout<<"StopFragmenting MinimalStringMass string->Mass() "<<MinimalStringMass<<" "<<string->Mass()<<G4endl; 
  G4cout<<"StopFragmenting - Yes/No "<<Result<<G4endl;
#endif	
		return Result;
	}
}

//-----------------------------------------------------------------------------
G4KineticTrack * G4LundStringFragmentation::Splitup(
		        G4FragmentingString *string, 
			G4FragmentingString *&newString)
{
#ifdef debug_LUNDfragmentation
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

// Uzhi 23 Dec. 2016
       G4double StringMass=string->Mass();
       G4double HadronMass=0.;
       G4double MinMass   =0.;
       G4int Iter(0), MaxIter(1);//000);   // Uzhi Jan. 2017
       do {
	    if (string->DecayIsQuark())
	    {
	     HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);
	    } else {
	     HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
	    }

#ifdef debug_LUNDfragmentation
G4cout<<"Iter "<<Iter<<G4endl;
  G4cout<<"The parton "<<string->GetDecayParton()->GetPDGEncoding()<<" "
        <<" produces hadron "<<HadronDefinition->GetParticleName()
        <<" and is transformed to "<<newStringEnd->GetPDGEncoding()<<G4endl;
  G4cout<<"The side of the string decay Left/Right (1/-1) "<<SideOfDecay<<G4endl;
#endif
// create new String from old, ie. keep Left and Right order, but replace decay

            if ( newString ) delete newString;  //AR-09Mar2017 Fixed memory leak

	    newString=new G4FragmentingString(*string,newStringEnd); // To store possible
                                                                     // quark containt of new string
	    SetMinimalStringMass(newString); MinMass = MinimalStringMass;
	    HadronMass = HadronDefinition->GetPDGMass();

	    if(Iter >= MaxIter) {G4KineticTrack * Hadron =0; return Hadron;}
	    Iter++;

//Uzhi Dec. 2016 G4cout<<"HadronMass + MinMass > StringMass "<<HadronMass<<" "<<MinMass<<" "<<HadronMass + MinMass<<" "<< StringMass<<G4endl;

       } while(HadronMass + MinMass > StringMass);

#ifdef debug_LUNDfragmentation
  G4cout<<"An attempt to determine its energy (SplitEandP)"<<G4endl;
#endif
       G4LorentzVector* HadronMomentum=SplitEandP(HadronDefinition, string, newString);

       delete newString; newString=0;
	
       G4KineticTrack * Hadron =0;
       if ( HadronMomentum != 0 ) {

#ifdef debug_LUNDfragmentation
  G4cout<<"The attempt was successful"<<G4endl;
#endif
	   G4ThreeVector   Pos;
	   Hadron = new G4KineticTrack(HadronDefinition, 0,Pos, *HadronMomentum);

           if ( newString ) delete newString;  //AR-09Mar2017 Fixed memory leak

	   newString=new G4FragmentingString(*string,newStringEnd,
	   				HadronMomentum);
	   delete HadronMomentum;
       }
       else
       {

#ifdef debug_LUNDfragmentation
  G4cout<<"The attempt was not successful !!!"<<G4endl;
#endif
       }

#ifdef debug_LUNDfragmentation
  G4cout<<"End SplitUP (G4VLongitudinalStringDecay) ====================="<<G4endl;
#endif

       return Hadron;
}

// Uzhi June 2014 Insert from G4ExcitedStringDecay.cc
//-----------------------------------------------------------------------------
G4ParticleDefinition * G4LundStringFragmentation::DiQuarkSplitup(
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
      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decayQuark);
      return had;
//      return hadronizer->Build(QuarkPair.first, decayQuark);

   } else {
   //... Diquark does not break

      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1;
                        // if we have a diquark, we need quark)

      G4double StrSup=GetStrangeSuppress();
      StrangeSuppress=0.43;   //0.42 0.38
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      StrangeSuppress=StrSup;

      created = QuarkPair.second;

      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
      return had;
   }
}

//-----------------------------------------------------------------------------
G4LorentzVector * G4LundStringFragmentation::SplitEandP(G4ParticleDefinition * pHadron,
                                                        G4FragmentingString *   string, 
                                                        G4FragmentingString * newString)
{ 
	G4LorentzVector String4Momentum=string->Get4Momentum();
	G4double StringMT2=string->MassT2();
	G4double StringMT =std::sqrt(StringMT2);

	G4double HadronMass = pHadron->GetPDGMass();

	SetMinimalStringMass(newString);

#ifdef debug_LUNDfragmentation
  G4cout<<G4endl<<"Start LUND SplitEandP "<<G4endl;
  G4cout<<"String 4 mom, String M and Mt "<<String4Momentum<<" "<<String4Momentum.mag()<<" "<<std::sqrt(StringMT2)<<G4endl;
  G4cout<<"Hadron "<<pHadron->GetParticleName()<<G4endl;
  G4cout<<"HadM MinimalStringMassLeft StringM hM+sM "<<HadronMass<<" "<<MinimalStringMass<<" "
        <<String4Momentum.mag()<<" "<<HadronMass+MinimalStringMass<<G4endl;
#endif

	if(HadronMass + MinimalStringMass > string->Mass()) 
	{
#ifdef debug_LUNDfragmentation
  G4cout<<"Mass of the string is not sufficient to produce the hadron!"<<G4endl;
#endif
	 return 0;
	}// have to start all over!

	String4Momentum.setPz(0.);
	G4ThreeVector StringPt=String4Momentum.vect();

	// calculate and assign hadron transverse momentum component HadronPx and HadronPy
	G4ThreeVector HadronPt    , RemSysPt; 
	G4double      HadronMassT2, ResidualMassT2;
	G4double HadronMt, Pt, Pt2, phi;             // Uzhi Nov. 2016 Introduction of Mt distribution
        //...  sample Pt of the hadron
        G4int attempt=0;
        do
        {
         attempt++; if(attempt > StringLoopInterrupt) {return 0;}  // Uzhi May 2016
                                                     // Uzhi Nov. 2016 Introduction of Mt distribution
	 HadronMt = HadronMass - 250.0*G4Log(G4UniformRand()); //200 -> 175
	 Pt2 = sqr(HadronMt)-sqr(HadronMass); Pt=std::sqrt(Pt2);
	 phi = 2.*pi*G4UniformRand();
	 G4ThreeVector SampleQuarkPtw= G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
                                                     // Uzhi Nov. 2016 Introduction of Mt distribution
         HadronPt =SampleQuarkPtw  + string->DecayPt();	// SampleQuarkPt()   // Uzhi Nov. 2016
         HadronPt.setZ(0);
         RemSysPt = StringPt - HadronPt;

         HadronMassT2 = sqr(HadronMass) + HadronPt.mag2();
         ResidualMassT2=sqr(MinimalStringMass) + RemSysPt.mag2();

        } while(std::sqrt(HadronMassT2) + std::sqrt(ResidualMassT2) > StringMT);

	//...  sample z to define hadron longitudinal momentum and energy
	//... but first check the available phase space

	G4double Pz2 = (sqr(StringMT2 - HadronMassT2 - ResidualMassT2) -
			4*HadronMassT2 * ResidualMassT2)/4./StringMT2;

	if(Pz2 < 0 ) {return 0;}          // have to start all over!

	//... then compute allowed z region  z_min <= z <= z_max

	G4double Pz = std::sqrt(Pz2);
	G4double zMin = (std::sqrt(HadronMassT2+Pz2) - Pz)/std::sqrt(StringMT2);
//	G4double zMin = (std::sqrt(HadronMassT2+Pz2) - 0.)/std::sqrt(StringMT2); // Uzhi 2014 For testing purposes
	G4double zMax = (std::sqrt(HadronMassT2+Pz2) + Pz)/std::sqrt(StringMT2);

	if (zMin >= zMax) return 0;		// have to start all over!

	G4double z = GetLightConeZ(zMin, zMax, 
			string->GetDecayParton()->GetPDGEncoding(), pHadron,
			HadronPt.x(), HadronPt.y());

	//... now compute hadron longitudinal momentum and energy
	// longitudinal hadron momentum component HadronPz

	HadronPt.setZ(0.5* string->GetDecayDirection() *
		      (z * string->LightConeDecay() - 
					HadronMassT2/(z * string->LightConeDecay())));
	G4double HadronE  = 0.5* (z * string->LightConeDecay() +
					HadronMassT2/(z * string->LightConeDecay()));

	G4LorentzVector * a4Momentum= new G4LorentzVector(HadronPt,HadronE);

#ifdef debug_LUNDfragmentation
 G4cout<<G4endl<<" string->GetDecayDirection() "<<string->GetDecayDirection()<<G4endl<<G4endl;
  G4cout<<"string->LightConeDecay() "<<string->LightConeDecay()<<G4endl;
  G4cout<<"HadronPt,HadronE "<<HadronPt<<" "<<HadronE<<G4endl;
  G4cout<<"String4Momentum "<<String4Momentum<<G4endl;
  G4cout<<"Out of LUND SplitEandP "<<G4endl<<G4endl;
#endif
        /* Adding to RestMomentum */
//        if((pHadron->GetPDGIsospin() >= 1.5 ) && (pHadron->GetBaryonNumber() >= 1.))
//	{RestMomentum +=G4LorentzVector(-a4Momentum->getX(),-a4Momentum->getY(),0.,0.);} // Uzhi June 2016
//
	return a4Momentum;
}

//-----------------------------------------------------------------------------------------
G4double G4LundStringFragmentation::GetLightConeZ(G4double zmin, G4double zmax, 
		G4int PDGEncodingOfDecayParton,
		G4ParticleDefinition* pHadron,
		G4double Px, G4double Py)
{
	G4double Mass = pHadron->GetPDGMass();
        G4int HadronEncoding=std::abs(pHadron->GetPDGEncoding());

	G4double Mt2 = Px*Px + Py*Py + Mass*Mass;

	G4double  Alund, Blund;
	G4double zOfMaxyf(0.), maxYf(1.), z(0.), yf(1.);

	if(!((std::abs(PDGEncodingOfDecayParton) > 1000) && (HadronEncoding > 1000)))
	{    // ---------------- Quark fragmentation  and qq-> meson ----------------------
Alund=1.;
          Blund=0.7/GeV/GeV;
//Mt          if((HadronEncoding == 321) || (HadronEncoding == 130) || (HadronEncoding == 310) ) Blund*=0.333;  // Uzhi Nov. 2016
//          if(std::abs(pHadron->GetBaryonNumber()) != 0) Alund=4.; // Blund/=8.; //Blund*=0.125;  // Uzhi Nov. 2016
	//    If alund get restored, you MUST adapt the calculation of zOfMaxyf.
	//    const G4double  Alund = 1;

G4double BMt2 = Blund*Mt2;
if(Alund == 1.0)
{	   zOfMaxyf=BMt2/(Blund*Mt2 + 1.);}
else
{  
   zOfMaxyf = ((1.0+BMt2) - std::sqrt(sqr(1.0-BMt2) + 4.0*BMt2*Alund))/2.0/(1.-Alund);
}

//G4cout<<"Alund Blund "<<Alund<<" "<<Blund<<G4endl;

if(zOfMaxyf < zmin) {zOfMaxyf=zmin;}               // Uzhi 8 Feb. 2016
if(zOfMaxyf > zmax) {zOfMaxyf=zmax;}
	   maxYf=(1-zOfMaxyf)/zOfMaxyf * G4Exp(-Blund*Mt2/zOfMaxyf);
//G4cout<<"zmin amax Zmax "<<zmin<<" "<<zmax<<" "<<zOfMaxyf<<G4endl;
//G4cout<<"BMt2 "<<BMt2; G4cout<<G4Exp(-BMt2/zmin)<<G4endl;
           const G4int maxNumberOfLoops = 1000;
           G4int loopCounter = 0;
	   do
	   {
		z = zmin + G4UniformRand()*(zmax-zmin);
//		yf = (1-z)/z * G4Exp(-Blund*Mt2/z);
		yf = G4Pow::GetInstance()->powA(1.0-z,Alund)/z*G4Exp(-BMt2/z);
	   }
	   while ( (G4UniformRand()*maxYf > yf) && ++loopCounter < maxNumberOfLoops );
           if ( loopCounter >= maxNumberOfLoops ) {
             z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
           }
	   return z;
        }

	if(std::abs(PDGEncodingOfDecayParton) > 1000)         // Uzhi May. 2016
	{
		G4double an=2.5; 
		if(pHadron->GetPDGIsospin() > 0.5) an=0.5;//0.75;
an +=(sqr(Px)+sqr(Py))/sqr(GeV)-0.5;
		z=zmin + (zmax-zmin)*G4Pow::GetInstance()->powA(G4UniformRand(),1./an);
	}

	return z;

}


// ####################################################

//----------------------------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::SplitLast(G4FragmentingString * string,
                                            G4KineticTrackVector * LeftVector,
                                            G4KineticTrackVector * RightVector)
{
	//... perform last cluster decay
// // Uzhi Dec. 2016
#ifdef debug_LUNDfragmentation
  G4cout<<G4endl<<"Split last-----------------------------------------"<<G4endl;
  G4cout<<"string->LightConePlus()  "<<string->LightConePlus() <<" "
        <<"string->LightConeMinus() "<<string->LightConeMinus()<<" "
        <<"string->LightConeDecay() "<<string->LightConeDecay()<<" "
        <<"string->GetDecayDirection() "<<string->GetDecayDirection()<<G4endl<<G4endl;
G4cout<<"String4mom "<<string->GetPstring()<<G4endl;
G4cout<<"Enter 4mom "<<string->Get4Momentum()<<G4endl;
#endif
// // Uzhi Dec.2016
	G4LorentzVector Str4Mom=string->Get4Momentum();
	G4ThreeVector ClusterVel=string->Get4Momentum().boostVector();
	G4double StringMass=string->Mass();

	G4ParticleDefinition * LeftHadron(0), * RightHadron(0);

        NumberOf_FS=0;
	for(G4int i=0; i<35; i++) {FS_Weight[i]=0.;}

	G4int sampledState = 0;

#ifdef debug_LUNDfragmentation
  G4cout<<"StrMass "<<StringMass<<" q's "
        <<string->GetLeftParton()->GetParticleName()<<" "
        <<string->GetRightParton()->GetParticleName()<<G4endl;
#endif

	string->SetLeftPartonStable(); // to query quark contents..

	if (string->FourQuarkString() )
	{
		// The string is qq-qqbar type. Diquarks are on the string ends
		//G4cout<<"The string is qq-qqbar type. Diquarks are on the string ends"<<G4endl;
	        if(StringMass-MinimalStringMass < 0.)
		{
			if (! Diquark_AntiDiquark_belowThreshold_lastSplitting(string, LeftHadron, RightHadron) ) 
                        {
				return false;
                        }
		} else
		{
			Diquark_AntiDiquark_aboveThreshold_lastSplitting(string, LeftHadron, RightHadron);

			if(NumberOf_FS == 0) return false;

                        sampledState = SampleState();
			if(string->GetLeftParton()->GetPDGEncoding() < 0)    // Uzhi June 2016
			{
				LeftHadron =FS_LeftHadron[sampledState];
				RightHadron=FS_RightHadron[sampledState];
			} else
			{
				LeftHadron =FS_RightHadron[sampledState];
				RightHadron=FS_LeftHadron[sampledState];
			}
		}
        } else
	{
		if (string->DecayIsQuark() && string->StableIsQuark() )
		{       //... there are quarks on cluster ends
#ifdef debug_LUNDfragmentation
  G4cout<<"Q Q string LastSplit"<<G4endl;
#endif
			Quark_AntiQuark_lastSplitting(string, LeftHadron, RightHadron);
		
			if(NumberOf_FS == 0) return false;
                	sampledState = SampleState();
			if(string->GetLeftParton()->GetPDGEncoding() < 0)            // Uzhi June 2016
			{
				LeftHadron =FS_RightHadron[sampledState];
				RightHadron=FS_LeftHadron[sampledState];
			} else
			{
				LeftHadron =FS_LeftHadron[sampledState];
				RightHadron=FS_RightHadron[sampledState];
			}
		} else 
		{       //... there is a Diquark on one of the cluster ends
#ifdef debug_LUNDfragmentation
  G4cout<<"DiQ Q string Last Split"<<G4endl;
#endif

			Quark_Diquark_lastSplitting(string, LeftHadron, RightHadron);

			if(NumberOf_FS == 0) return false;                           
                	sampledState = SampleState();

			if(string->GetLeftParton()->GetParticleSubType() == "quark") // Uzhi June 2016
			{
				LeftHadron =FS_LeftHadron[sampledState];
				RightHadron=FS_RightHadron[sampledState];
			} else 
			{
				LeftHadron =FS_RightHadron[sampledState];
				RightHadron=FS_LeftHadron[sampledState];
			}			
		}

/*		
		if(NumberOf_FS == 0) return false;                                   // Uzhi June 2016
                G4int sampledState = SampleState();
		LeftHadron =FS_LeftHadron[sampledState];
		RightHadron=FS_RightHadron[sampledState];
*/
#ifdef debug_LUNDfragmentation
  G4cout<<"Selected LeftHad RightHad "<<sampledState<<" "
        <<LeftHadron->GetParticleName()<<" "<<RightHadron->GetParticleName()<<G4endl;
#endif

	}  // End of if(!string->FourQuarkString())

	G4LorentzVector  LeftMom, RightMom;
	G4ThreeVector    Pos;

	Sample4Momentum(&LeftMom,  LeftHadron->GetPDGMass(),
			&RightMom, RightHadron->GetPDGMass(),
			StringMass);

	LeftMom.boost(ClusterVel);
	RightMom.boost(ClusterVel);
/*
G4cout<<G4endl<<"StringMass "<<StringMass<<G4endl;
G4cout<<LeftMom<<G4endl;
G4cout<<RightMom<<G4endl;
{G4int Uzhi; G4cin>>Uzhi;}
*/
	LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
	RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

	return true;

}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Diquark_AntiDiquark_belowThreshold_lastSplitting(G4FragmentingString * & string,
                                                 G4ParticleDefinition * & LeftHadron,
                                                 G4ParticleDefinition * & RightHadron)
{
	G4double StringMass   = string->Mass();

	G4int cClusterInterrupt = 0;
	do
	{
		if (cClusterInterrupt++ >= ClusterLoopInterrupt)
		{
			return false;
		}

		G4int LeftQuark1= string->GetLeftParton()->GetPDGEncoding()/1000;
		G4int LeftQuark2=(string->GetLeftParton()->GetPDGEncoding()/100)%10;

		G4int RightQuark1= string->GetRightParton()->GetPDGEncoding()/1000;
		G4int RightQuark2=(string->GetRightParton()->GetPDGEncoding()/100)%10;

		if(G4UniformRand()<0.5)
		{
			LeftHadron =hadronizer->Build(FindParticle( LeftQuark1),
						      FindParticle(RightQuark1));
			RightHadron=hadronizer->Build(FindParticle( LeftQuark2),
						      FindParticle(RightQuark2));
		} else
		{
			LeftHadron =hadronizer->Build(FindParticle( LeftQuark1),
						      FindParticle(RightQuark2));
			RightHadron=hadronizer->Build(FindParticle( LeftQuark2),
						      FindParticle(RightQuark1));
		}

		//... repeat procedure, if mass of cluster is too low to produce hadrons
		//... ClusterMassCut = 0.15*GeV model parameter
	}
	while ((StringMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()));

  	return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Diquark_AntiDiquark_aboveThreshold_lastSplitting(G4FragmentingString * & string,
                                                 G4ParticleDefinition * & LeftHadron,
                                                 G4ParticleDefinition * & RightHadron)
{
	// StringMass-MinimalStringMass > 0. Creation of 2 baryons is possible ----

	G4double StringMass   = string->Mass();
	G4double StringMassSqr= sqr(StringMass); 
	G4ParticleDefinition * Di_Quark;
	G4ParticleDefinition * Anti_Di_Quark;

	if(string->GetLeftParton()->GetPDGEncoding() < 0)
	{
		Anti_Di_Quark   =string->GetLeftParton();
		Di_Quark=string->GetRightParton();
	} else
	{
		Anti_Di_Quark   =string->GetRightParton();
		Di_Quark=string->GetLeftParton();
	}

	G4int IDAnti_di_quark    =Anti_Di_Quark->GetPDGEncoding();
	G4int AbsIDAnti_di_quark =std::abs(IDAnti_di_quark);
	G4int IDdi_quark         =Di_Quark->GetPDGEncoding();
	G4int AbsIDdi_quark      =std::abs(IDdi_quark);

	G4int ADi_q1=AbsIDAnti_di_quark/1000;
	G4int ADi_q2=(AbsIDAnti_di_quark-ADi_q1*1000)/100;

	G4int Di_q1=AbsIDdi_quark/1000;
	G4int Di_q2=(AbsIDdi_quark-Di_q1*1000)/100;

	NumberOf_FS=0;
	for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
	{
		G4int StateADiQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do  // while(Meson[AbsIDquark-1][ProdQ-1][StateQ]<>0);
		{
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(
							-Baryon[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]);
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			//G4cout<<"Anti Bar "<<LeftHadron->GetParticleName()<<G4endl;

			G4int StateDiQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do // while(Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]<>0);
			{
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(
								+Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if(StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 34 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Diquark_AntiDiquark_aboveThreshold_lastSplitting ",
                                                       "HAD_LUND_001", JustWarning, ed );
                                          NumberOf_FS = 34;
                                        }

					G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),
								sqr(RightHadronMass));
					//FS_Psqr=1.;
					FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*FS_Psqr*
							       BaryonWeight[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]*
							       BaryonWeight[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]*
							       Prob_QQbar[ProdQ-1];

					FS_LeftHadron[NumberOf_FS] = LeftHadron;
					FS_RightHadron[NumberOf_FS]= RightHadron;

					NumberOf_FS++;
				} // End of if(StringMass > LeftHadronMass + RightHadronMass)

				StateDiQ++;

			} while( (Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]!=0) && 
                                 ++internalLoopCounter < maxNumberOfInternalLoops );
                        if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
                          return false;
                        }
 
			StateADiQ++;
		} while( (Baryon[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]!=0) &&
                         ++loopCounter < maxNumberOfLoops ); 
                if ( loopCounter >= maxNumberOfLoops ) {
                  return false;
                }
	} // End of for(G4int ProdQ=1; ProdQ < 4; ProdQ++)

  	return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Quark_Diquark_lastSplitting(G4FragmentingString * & string,
                            G4ParticleDefinition * & LeftHadron,
                            G4ParticleDefinition * & RightHadron)
{
	G4double StringMass   = string->Mass();
	G4double StringMassSqr= sqr(StringMass);

	G4ParticleDefinition * Di_Quark;
	G4ParticleDefinition * Quark;

	if(string->GetLeftParton()->GetParticleSubType()== "quark")
	{
		Quark   =string->GetLeftParton();
		Di_Quark=string->GetRightParton();
	} else
	{
		Quark   =string->GetRightParton();
		Di_Quark=string->GetLeftParton();
	}

	G4int IDquark        =Quark->GetPDGEncoding();
	G4int AbsIDquark     =std::abs(IDquark);
	G4int IDdi_quark   =Di_Quark->GetPDGEncoding();
	G4int AbsIDdi_quark=std::abs(IDdi_quark);
	G4int Di_q1=AbsIDdi_quark/1000;
	G4int Di_q2=(AbsIDdi_quark-Di_q1*1000)/100;

	G4int              SignDiQ= 1;
	if(IDdi_quark < 0) SignDiQ=-1;

	NumberOf_FS=0;
	for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
	{
		G4int SignQ;
		if(IDquark > 0)
		{                                   SignQ=-1;
			if(IDquark == 2)                   SignQ= 1;
			if((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0
			if((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar
		} else
		{
			SignQ= 1;
			if(IDquark == -2)                  SignQ=-1;
			if((IDquark ==-1) && (ProdQ == 3)) SignQ=-1; // K0bar
			if((IDquark ==-3) && (ProdQ == 1)) SignQ= 1; // K0
		}

		if(AbsIDquark == ProdQ)            SignQ= 1;

		G4int StateQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do  // while(Meson[AbsIDquark-1][ProdQ-1][StateQ]<>0);
		{
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignQ*
							Meson[AbsIDquark-1][ProdQ-1][StateQ]);
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			G4int StateDiQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do // while(Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]<>0);
			{
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignDiQ*
								Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if(StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 34 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Quark_Diquark_lastSplitting ",
                                                       "HAD_LUND_002", JustWarning, ed );
                                          NumberOf_FS = 34;
                                        }

					G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),
								sqr(RightHadronMass));
					FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*
							       MesonWeight[AbsIDquark-1][ProdQ-1][StateQ]*
							       BaryonWeight[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]*
							       Prob_QQbar[ProdQ-1];

					FS_LeftHadron[NumberOf_FS] = LeftHadron;
					FS_RightHadron[NumberOf_FS]= RightHadron;

					NumberOf_FS++;
				} // End of if(StringMass > LeftHadronMass + RightHadronMass)

				StateDiQ++;

			} while( (Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]!=0) &&
                                 ++internalLoopCounter < maxNumberOfInternalLoops );
                        if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
                          return false;
                        }

			StateQ++;
		} while( (Meson[AbsIDquark-1][ProdQ-1][StateQ]!=0) &&
                         ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */

                  if ( loopCounter >= maxNumberOfLoops ) {
                    return false;
                  }
	} // End of for(G4int ProdQ=1; ProdQ < 4; ProdQ++)

	return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Quark_AntiQuark_lastSplitting(G4FragmentingString * & string,
                              G4ParticleDefinition * & LeftHadron,
                              G4ParticleDefinition * & RightHadron)
{
	G4double StringMass   = string->Mass();
	G4double StringMassSqr= sqr(StringMass);

	G4ParticleDefinition * Quark;
	G4ParticleDefinition * Anti_Quark;

	if(string->GetLeftParton()->GetPDGEncoding()>0)
	{
		Quark     =string->GetLeftParton();
		Anti_Quark=string->GetRightParton();
	} else
	{
		Quark     =string->GetRightParton();
		Anti_Quark=string->GetLeftParton();
	}

	G4int IDquark        =Quark->GetPDGEncoding();
	G4int AbsIDquark     =std::abs(IDquark);
	G4int IDanti_quark   =Anti_Quark->GetPDGEncoding();
	G4int AbsIDanti_quark=std::abs(IDanti_quark);

	NumberOf_FS=0;
	for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
	{
		G4int                              SignQ=-1;
		if(IDquark == 2)                   SignQ= 1;
		if((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0
		if((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar
		if(IDquark == ProdQ)               SignQ= 1;

		G4int                                   SignAQ= 1;
		if(IDanti_quark == -2)                  SignAQ=-1;
		if((IDanti_quark ==-1) && (ProdQ == 3)) SignAQ=-1; // K0bar
		if((IDanti_quark ==-3) && (ProdQ == 1)) SignAQ= 1; // K0
		if(AbsIDanti_quark == ProdQ)            SignAQ= 1;

		G4int StateQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do  // while(Meson[AbsIDquark-1][ProdQ-1][StateQ]<>0);
		{
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignQ*
						       Meson[AbsIDquark-1][ProdQ-1][StateQ]);
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			G4int StateAQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do // while(Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]<>0);
			{
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignAQ*
								Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]);
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if(StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 34 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Quark_AntiQuark_lastSplitting ",
                                                       "HAD_LUND_003", JustWarning, ed );
                                          NumberOf_FS = 34;
                                        }

                                        G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),
								sqr(RightHadronMass));
					//FS_Psqr=1.;
					FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*
							       MesonWeight[AbsIDquark-1][ProdQ-1][StateQ]*
							       MesonWeight[AbsIDanti_quark-1][ProdQ-1][StateAQ]*
							       Prob_QQbar[ProdQ-1];

					if(string->GetLeftParton()->GetPDGEncoding()>0)
					{
						FS_LeftHadron[NumberOf_FS] = RightHadron;
						FS_RightHadron[NumberOf_FS]= LeftHadron;
					} else
					{
						FS_LeftHadron[NumberOf_FS] = LeftHadron;
						FS_RightHadron[NumberOf_FS]= RightHadron;
					}

					NumberOf_FS++;

				} // End of if(StringMass > LeftHadronMass + RightHadronMass)

				StateAQ++;
			} while( (Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]!=0) &&
                                 ++internalLoopCounter < maxNumberOfInternalLoops );
                          if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
                            return false;
                          }

			StateQ++;
		} while( (Meson[AbsIDquark-1][ProdQ-1][StateQ]!=0) &&
                         ++loopCounter < maxNumberOfLoops );
                  if ( loopCounter >= maxNumberOfLoops ) {
                    return false;
                  }
	} // End of for(G4int ProdQ=1; ProdQ < 4; ProdQ++)

	return true;
}

//----------------------------------------------------------------------------------------------------------
G4int G4LundStringFragmentation::SampleState(void) 
{
        if ( NumberOf_FS > 34 ) {
          G4ExceptionDescription ed;
          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
          G4Exception( "G4LundStringFragmentation::SampleState ", "HAD_LUND_004", JustWarning, ed );
          NumberOf_FS = 34;
        }

	G4double SumWeights=0.;

	for(G4int i=0; i<NumberOf_FS; i++) {SumWeights+=FS_Weight[i];}// G4cout<<i<<" "<<FS_Weight[i]<<G4endl;}

	G4double ksi=G4UniformRand();
	G4double Sum=0.;
	G4int indexPosition = 0;

	for(G4int i=0; i<NumberOf_FS; i++)
	{
		Sum+=(FS_Weight[i]/SumWeights);
		indexPosition=i;
		if(Sum >= ksi) break;
	}
	return indexPosition;
}

//----------------------------------------------------------------------------------------------------------
void G4LundStringFragmentation::Sample4Momentum(G4LorentzVector*     Mom, G4double     Mass, 
                                                G4LorentzVector* AntiMom, G4double AntiMass, 
                                                G4double InitialMass) 
{
	// ------ Sampling of momenta of 2 last produced hadrons --------------------
	G4ThreeVector Pt;
	G4double MassMt, AntiMassMt;
	G4double AvailablePz, AvailablePz2;
        G4double ProbIsotropy = sqr(sqr(938.0/InitialMass));                     // Uzhi May 2015
#ifdef debug_LUNDfragmentation
  G4cout<<"Sampling of momenta of 2 last produced hadrons ----------------"<<G4endl;
  G4cout<<"Masses "<<InitialMass<<" "<<Mass<<" "<<AntiMass<<" ProbIsotropy "<<ProbIsotropy<<G4endl;
#endif

	G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) -    // Uzhi June 2016
 			sqr(2.*Mass*AntiMass);
	G4double Pabs = (r_val > 0.)? std::sqrt(r_val)/(2.*InitialMass) : 0;

	if((Mass > 930. || AntiMass > 930.) && (G4UniformRand() < ProbIsotropy))
	{       //If there is a baryon
		// ----------------- Isotropic decay ------------------------------------
		//...              sample unit vector
		G4double pz =1. - 2.*G4UniformRand();
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
	else
	{
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do
		{
			// GF 22-May-09, limit sampled pt to allowed range
/* 												// Uzhi June 2016
			G4double termD = InitialMass*InitialMass -Mass*Mass - AntiMass*AntiMass;
			G4double termab = 4*sqr(Mass*AntiMass);
			G4double termN = 2*termD + 4*Mass*Mass + 4*AntiMass*AntiMass;
			G4double pt2max=(termD*termD - termab )/ termN ;
			//G4cout<<"Anis "<<pt2max<<" "<<(termD*termD-termab)/(4.*InitialMass*InitialMass)<<G4endl;

			Pt=SampleQuarkPt(std::sqrt(pt2max)); Pt.setZ(0); G4double Pt2=Pt.mag2();
*/
			Pt=SampleQuarkPt(Pabs); Pt.setZ(0); G4double Pt2=Pt.mag2();

			MassMt    = std::sqrt(    Mass *     Mass + Pt2);
			AntiMassMt= std::sqrt(AntiMass * AntiMass + Pt2);

//			AvailablePz2= sqr(InitialMass*InitialMass - MassMt2 - AntiMassMt2) -
//					4.*MassMt2*AntiMassMt2;
//		}
//		while( (AvailablePz2 < 0.) &&  // GF will occur only for numerical precision problem with limit in sampled pt
//                    ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */ // Uzhi June 2016

		}
		while( (InitialMass < MassMt + AntiMassMt) && ++loopCounter < maxNumberOfLoops ); 

                if ( loopCounter >= maxNumberOfLoops ) {
                  AvailablePz2 = 0.0;
                }

		AvailablePz2= sqr(InitialMass*InitialMass - sqr(MassMt) - sqr(AntiMassMt)) -
				4.*sqr(MassMt*AntiMassMt);

		AvailablePz2 /=(4.*InitialMass*InitialMass);
		AvailablePz = std::sqrt(AvailablePz2);

		G4double Px=Pt.getX();
		G4double Py=Pt.getY();

		Mom->setPx(Px); Mom->setPy(Py); Mom->setPz(AvailablePz);
		Mom->setE(std::sqrt(sqr(MassMt)+AvailablePz2));

		AntiMom->setPx(-Px); AntiMom->setPy(-Py); AntiMom->setPz(-AvailablePz);
		AntiMom->setE (std::sqrt(sqr(AntiMassMt)+AvailablePz2));
	}
}

//------------------------------------------------------------------------
G4double G4LundStringFragmentation::lambda(G4double S, G4double m1_Sqr, G4double m2_Sqr)
{ 
	G4double lam = sqr(S - m1_Sqr - m2_Sqr) - 4.*m1_Sqr*m2_Sqr;
	return lam;
}

//-----------------------------------------------------------------------
void G4LundStringFragmentation::SetMinMasses()
{
// ------ For estimation of a minimal string mass ---------------
    Mass_of_light_quark    =140.*MeV;
    Mass_of_heavy_quark    =500.*MeV;
    Mass_of_string_junction=720.*MeV;

    G4double minMQQbarStr[3][3] ={ {350.0, 350.0, 710.0},  //DDbar, DUbar, DSbar in MeV                 // Uzhi July 2016
                                   {350.0, 350.0, 710.0},  //UDbar, UUbar, USbar in Mev
                                   {710.0, 710.0,1070.0 }};//SDbar, SUbar, SSbar in MeV
    for(G4int i=0; i<3; i++){ for(G4int j=0; j<3; j++){minMassQQbarStr[i][j]=minMQQbarStr[i][j];};}; // Uzhi July 2016

 
   G4double minMQDiQStr[3][3][3] = {{{1160., 1160., 1340.}, {1160., 1160., 1340.}, {1340., 1340., 1540.},}, //d-dd, d-du, d-ds, d-ud, d-uu, d-us, d-sd, d-su, d-ss // Uzhi July 2016
                                    {{1160., 1160., 1340.}, {1160., 1160., 1340.}, {1340., 1340., 1540.},}, //u-dd, u-du, u-ds, u-ud, u-uu, u-us, u-sd, u-su, u-ss
                                    {{1520., 1520., 1690.}, {1520., 1520., 1690.}, {1690., 1690., 1890. }}};//s-dd, s-du, s-ds, s-ud, s-uu, s-us, s-sd, s-su, s-ss
   for(G4int i=0; i<3; i++){ for(G4int j=0; j<3; j++){ for(G4int k=0; k<3; k++){minMassQDiQStr[i][j][k]=minMQDiQStr[i][j][k];};};};                      // Uzhi July 2016

// ------ An estimated minimal string mass ----------------------
    MinimalStringMass  = 0.;
    MinimalStringMass2 = 0.;

// For treating of small string decays
   for(G4int i=0; i<3; i++)
   {  for(G4int j=0; j<3; j++)
      {  for(G4int k=0; k<6; k++)
         {
           Meson[i][j][k]=0; MesonWeight[i][j][k]=0.;
         }
      }
   }
//--------------------------
         Meson[0][0][0]=111;                       // dbar-d Pi0
   MesonWeight[0][0][0]=(1.-pspin_meson)*(1.-scalarMesonMix[0]);

         Meson[0][0][1]=221;                       // dbar-d Eta
   MesonWeight[0][0][1]=(1.-pspin_meson)*(scalarMesonMix[0]-scalarMesonMix[1]);

         Meson[0][0][2]=331;                       // dbar-d EtaPrime
   MesonWeight[0][0][2]=(1.-pspin_meson)*(scalarMesonMix[1]);

         Meson[0][0][3]=113;                       // dbar-d Rho0
   MesonWeight[0][0][3]=pspin_meson*(1.-vectorMesonMix[0]);

         Meson[0][0][4]=223;                       // dbar-d Omega
   MesonWeight[0][0][4]=pspin_meson*(vectorMesonMix[0]);
//--------------------------

         Meson[0][1][0]=211;                       // dbar-u Pi+
   MesonWeight[0][1][0]=(1.-pspin_meson);

         Meson[0][1][1]=213;                       // dbar-u Rho+
   MesonWeight[0][1][1]=pspin_meson;
//--------------------------

         Meson[0][2][0]=311;                      // dbar-s K0bar
   MesonWeight[0][2][0]=(1.-pspin_meson);

         Meson[0][2][1]=313;                       // dbar-s K*0bar
   MesonWeight[0][2][1]=pspin_meson;
//--------------------------
//--------------------------
         Meson[1][0][0]=211;                       // ubar-d Pi-
   MesonWeight[1][0][0]=(1.-pspin_meson);

         Meson[1][0][1]=213;                       // ubar-d Rho-
   MesonWeight[1][0][1]=pspin_meson;
//--------------------------

         Meson[1][1][0]=111;                       // ubar-u Pi0
   MesonWeight[1][1][0]=(1.-pspin_meson)*(1.-scalarMesonMix[0]);

         Meson[1][1][1]=221;                       // ubar-u Eta
   MesonWeight[1][1][1]=(1.-pspin_meson)*(scalarMesonMix[0]-scalarMesonMix[1]);

         Meson[1][1][2]=331;                       // ubar-u EtaPrime
   MesonWeight[1][1][2]=(1.-pspin_meson)*(scalarMesonMix[1]);

         Meson[1][1][3]=113;                       // ubar-u Rho0
   MesonWeight[1][1][3]=pspin_meson*(1.-vectorMesonMix[0]);

         Meson[1][1][4]=223;                       // ubar-u Omega
   //MesonWeight[1][1][4]=pspin_meson*(scalarMesonMix[0]);
   MesonWeight[1][1][4]=pspin_meson*(vectorMesonMix[0]);  // Uzhi 2015 scalar -> vector
//--------------------------

         Meson[1][2][0]=321;                      // ubar-s K-
   MesonWeight[1][2][0]=(1.-pspin_meson);

         Meson[1][2][1]=323;                      // ubar-s K*-bar -
   MesonWeight[1][2][1]=pspin_meson;
//--------------------------
//--------------------------

         Meson[2][0][0]=311;                       // sbar-d K0
   MesonWeight[2][0][0]=(1.-pspin_meson);

         Meson[2][0][1]=313;                       // sbar-d K*0
   MesonWeight[2][0][1]=pspin_meson;
//--------------------------

         Meson[2][1][0]=321;                        // sbar-u K+
   MesonWeight[2][1][0]=(1.-pspin_meson);

         Meson[2][1][1]=323;                       // sbar-u K*+
   MesonWeight[2][1][1]=pspin_meson;
//--------------------------

         Meson[2][2][0]=221;                       // sbar-s Eta
   MesonWeight[2][2][0]=(1.-pspin_meson)*(1.-scalarMesonMix[5]);

         Meson[2][2][1]=331;                       // sbar-s EtaPrime
   MesonWeight[2][2][1]=(1.-pspin_meson)*(1.-scalarMesonMix[5]);

         Meson[2][2][3]=333;                       // sbar-s EtaPrime
   MesonWeight[2][2][3]=pspin_meson*(vectorMesonMix[5]);
//--------------------------

   for(G4int i=0; i<3; i++)
   {  for(G4int j=0; j<3; j++)
      {  for(G4int k=0; k<3; k++)
         {  for(G4int l=0; l<4; l++)
            { Baryon[i][j][k][l]=0; BaryonWeight[i][j][k][l]=0.;}
         }
      }
   }

   G4double pspin_barion_in=pspin_barion;
   //pspin_barion=0.75;
//---------------------------------------
         Baryon[0][0][0][0]=1114;         // Delta-
   BaryonWeight[0][0][0][0]=1.;

//---------------------------------------
         Baryon[0][0][1][0]=2112;         // neutron
   BaryonWeight[0][0][1][0]=1.-pspin_barion;

         Baryon[0][0][1][1]=2114;         // Delta0
   BaryonWeight[0][0][1][1]=pspin_barion;

//---------------------------------------
         Baryon[0][0][2][0]=3112;         // Sigma-
   BaryonWeight[0][0][2][0]=1.-pspin_barion;

         Baryon[0][0][2][1]=3114;         // Sigma*-
   BaryonWeight[0][0][2][1]=pspin_barion;

//---------------------------------------
         Baryon[0][1][0][0]=2112;         // neutron
   BaryonWeight[0][1][0][0]=1.-pspin_barion;

         Baryon[0][1][0][1]=2114;         // Delta0
   BaryonWeight[0][1][0][1]=pspin_barion;

//---------------------------------------
         Baryon[0][1][1][0]=2212;         // proton
   BaryonWeight[0][1][1][0]=1.-pspin_barion;

         Baryon[0][1][1][1]=2214;         // Delta+
   BaryonWeight[0][1][1][1]=pspin_barion;

//---------------------------------------
         Baryon[0][1][2][0]=3122;         // Lambda
   BaryonWeight[0][1][2][0]=(1.-pspin_barion)*0.5;

         Baryon[0][1][2][1]=3212;         // Sigma0
   BaryonWeight[0][1][2][1]=(1.-pspin_barion)*0.5;

         Baryon[0][1][2][2]=3214;         // Sigma*0
   BaryonWeight[0][1][2][2]=pspin_barion;

//---------------------------------------
         Baryon[0][2][0][0]=3112;         // Sigma-
   BaryonWeight[0][2][0][0]=1.-pspin_barion;

         Baryon[0][2][0][1]=3114;         // Sigma*-
   BaryonWeight[0][2][0][1]=pspin_barion;

//---------------------------------------
         Baryon[0][2][1][0]=3122;         // Lambda
   BaryonWeight[0][2][1][0]=(1.-pspin_barion)*0.5;

         Baryon[0][2][1][1]=3212;         // Sigma0
   BaryonWeight[0][2][1][1]=(1.-pspin_barion)*0.5;

         Baryon[0][2][1][2]=3214;         // Sigma*0
   BaryonWeight[0][2][1][2]=pspin_barion;

//---------------------------------------
         Baryon[0][2][2][0]=3312;         // Theta-
   BaryonWeight[0][2][2][0]=1.-pspin_barion;

         Baryon[0][2][2][1]=3314;         // Theta*-
   BaryonWeight[0][2][2][1]=pspin_barion;

//---------------------------------------
//---------------------------------------
         Baryon[1][0][0][0]=2112;         // neutron
   BaryonWeight[1][0][0][0]=1.-pspin_barion;

         Baryon[1][0][0][1]=2114;         // Delta0
   BaryonWeight[1][0][0][1]=pspin_barion;

//---------------------------------------
         Baryon[1][0][1][0]=2212;         // proton
   BaryonWeight[1][0][1][0]=1.-pspin_barion;          

         Baryon[1][0][1][1]=2214;         // Delta+
   BaryonWeight[1][0][1][1]=pspin_barion;

//---------------------------------------
         Baryon[1][0][2][0]=3122;         // Lambda
   BaryonWeight[1][0][2][0]=(1.-pspin_barion)*0.5;

         Baryon[1][0][2][1]=3212;         // Sigma0
   BaryonWeight[1][0][2][1]=(1.-pspin_barion)*0.5;

         Baryon[1][0][2][2]=3214;         // Sigma*0
   BaryonWeight[1][0][2][2]=pspin_barion;

//---------------------------------------
         Baryon[1][1][0][0]=2212;         // proton
   BaryonWeight[1][1][0][0]=1.-pspin_barion;

         Baryon[1][1][0][1]=2214;         // Delta+
   BaryonWeight[1][1][0][1]=pspin_barion;

//---------------------------------------
         Baryon[1][1][1][0]=2224;         // Delta++
   BaryonWeight[1][1][1][0]=1.;

//---------------------------------------
         Baryon[1][1][2][0]=3222;         // Sigma+
   BaryonWeight[1][1][2][0]=1.-pspin_barion;

         Baryon[1][1][2][1]=3224;         // Sigma*+
   BaryonWeight[1][1][2][1]=pspin_barion;

//---------------------------------------
         Baryon[1][2][0][0]=3122;         // Lambda
   BaryonWeight[1][2][0][0]=(1.-pspin_barion)*0.5;

         Baryon[1][2][0][1]=3212;         // Sigma0
   BaryonWeight[1][2][0][1]=(1.-pspin_barion)*0.5;

         Baryon[1][2][0][2]=3214;         // Sigma*0
   BaryonWeight[1][2][0][2]=pspin_barion;

//---------------------------------------
         Baryon[1][2][1][0]=3222;         // Sigma+
   BaryonWeight[1][2][1][0]=1.-pspin_barion;

         Baryon[1][2][1][1]=3224;         // Sigma*+
   BaryonWeight[1][2][1][1]=pspin_barion;

//---------------------------------------
         Baryon[1][2][2][0]=3322;         // Theta0
   BaryonWeight[1][2][2][0]=1.-pspin_barion;

         Baryon[1][2][2][1]=3324;         // Theta*0
   BaryonWeight[1][2][2][1]=pspin_barion;

//---------------------------------------
//---------------------------------------
         Baryon[2][0][0][0]=3112;         // Sigma-
   BaryonWeight[2][0][0][0]=1.-pspin_barion;

         Baryon[2][0][0][1]=3114;         // Sigma*-
   BaryonWeight[2][0][0][1]=pspin_barion;

//---------------------------------------
         Baryon[2][0][1][0]=3122;         // Lambda
   BaryonWeight[2][0][1][0]=(1.-pspin_barion)*0.5;          

         Baryon[2][0][1][1]=3212;         // Sigma0
   BaryonWeight[2][0][1][1]=(1.-pspin_barion)*0.5; 

         Baryon[2][0][1][2]=3214;         // Sigma*0
   BaryonWeight[2][0][1][2]=pspin_barion;

//---------------------------------------
         Baryon[2][0][2][0]=3312;         // Sigma-
   BaryonWeight[2][0][2][0]=1.-pspin_barion;

         Baryon[2][0][2][1]=3314;         // Sigma*-
   BaryonWeight[2][0][2][1]=pspin_barion;

//---------------------------------------
         Baryon[2][1][0][0]=3122;         // Lambda
   BaryonWeight[2][1][0][0]=(1.-pspin_barion)*0.5;

         Baryon[2][1][0][1]=3212;         // Sigma0
   BaryonWeight[2][1][0][1]=(1.-pspin_barion)*0.5;

         Baryon[2][1][0][2]=3214;         // Sigma*0
   BaryonWeight[2][1][0][2]=pspin_barion;

//---------------------------------------
         Baryon[2][1][1][0]=3222;         // Sigma+
   BaryonWeight[2][1][1][0]=1.-pspin_barion;

         Baryon[2][1][1][1]=3224;         // Sigma*+
   BaryonWeight[2][1][1][1]=pspin_barion;

//---------------------------------------
         Baryon[2][1][2][0]=3322;         // Theta0
   BaryonWeight[2][1][2][0]=1.-pspin_barion;

         Baryon[2][1][2][1]=3324;         // Theta*0
   BaryonWeight[2][1][2][1]=pspin_barion;

//---------------------------------------
         Baryon[2][2][0][0]=3312;         // Theta-
   BaryonWeight[2][2][0][0]=1.-pspin_barion;

         Baryon[2][2][0][1]=3314;         // Theta*-
   BaryonWeight[2][2][0][1]=pspin_barion;

//---------------------------------------
         Baryon[2][2][1][0]=3322;         // Theta0
   BaryonWeight[2][2][1][0]=1.-pspin_barion;

         Baryon[2][2][1][1]=3324;         // Theta*0
   BaryonWeight[2][2][1][1]=pspin_barion;

//---------------------------------------
         Baryon[2][2][2][0]=3334;         // Omega
   BaryonWeight[2][2][2][0]=1.;

//---------------------------------------
   pspin_barion=pspin_barion_in;
   /*
	   for(G4int i=0; i<3; i++)
	   {  for(G4int j=0; j<3; j++)
		  {  for(G4int k=0; k<3; k++)
			 {  for(G4int l=0; l<4; l++)
				{ G4cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Baryon[i][j][k][l]<<G4endl;}
			 }
		  }
	   }
		G4int Uzhi;
		G4cin>>Uzhi;
    */

   G4double StrSup=GetStrangeSuppress();
   SetStrangenessSuppression(0.375);
   Prob_QQbar[0]=StrangeSuppress;         // Probability of ddbar production
   Prob_QQbar[1]=StrangeSuppress;         // Probability of uubar production
   Prob_QQbar[2]=1.0-2.*StrangeSuppress;  // Probability of ssbar production 
   SetStrangenessSuppression(StrSup);     // Uzhi June 2016 0.45->0.4573

//Uzhi Feb. 2016   SetStrangenessSuppression(0.45);       // Uzhi Nov. 2016
   for ( G4int i=0 ; i<35 ; i++ ) { 
     FS_LeftHadron[i] = 0;
     FS_RightHadron[i] = 0;
     FS_Weight[i] = 0.0; 
   }
//Uzhi Feb. 2016   SetStrangenessSuppression(StrSup);
   NumberOf_FS = 0;

}

// --------------------------------------------------------------
G4LundStringFragmentation::~G4LundStringFragmentation()
{}

//--------------------------------------------------------------------------------------
void G4LundStringFragmentation::SetMinimalStringMass(const G4FragmentingString  * const string)  
{
	G4double EstimatedMass=0.;

        G4int Qleft =std::abs(string->GetLeftParton()->GetPDGEncoding());   // Uzhi July 2016 Start
        G4int Qright=std::abs(string->GetRightParton()->GetPDGEncoding());

        if((Qleft < 4) && (Qright < 4)) {   // Q-Qbar string
          EstimatedMass=minMassQQbarStr[Qleft-1][Qright-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if((Qleft < 4) && (Qright > 1000)) {   // Q - DiQ string
          G4int q1=Qright/1000;
          G4int q2=(Qright/100)%10;
          EstimatedMass=minMassQDiQStr[Qleft-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if((Qleft > 1000) && (Qright < 4)) {   // DiQ - Q string
          G4int q1=Qleft/1000;
          G4int q2=(Qleft/100)%10;
          EstimatedMass=minMassQDiQStr[Qright-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }                                                                // Uzhi July 2016 End

//	DiQuark - Anti DiQuark string -----------------
	G4int Number_of_quarks=0;
        G4int Number_of_squarks=0;
        
	G4double StringM=string->Get4Momentum().mag();

#ifdef debug_LUNDfragmentation
//  G4cout<<"MinStringMass// Input String mass "<<string->Get4Momentum().mag()<<" Qleft "<<Qleft<<G4endl;
#endif

	if( Qleft > 1000)
	{
		Number_of_quarks+=2;
		G4int q1=Qleft/1000;
		if( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

		G4int q2=(Qleft/100)%10;
		if( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
//		EstimatedMass +=Mass_of_string_junction;
	}

#ifdef debug_LUNDfragmentation
//  G4cout<<"Min mass with Qleft "<<Qleft<<" "<<EstimatedMass<<G4endl;
#endif
//	G4int Qright=std::abs(string->GetRightParton()->GetPDGEncoding());

	if( Qright > 1000)
	{
		Number_of_quarks+=2;
		G4int q1=Qright/1000;
		if( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

		G4int q2=(Qright/100)%10;
		if( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
		//EstimatedMass +=Mass_of_string_junction;
	}

#ifdef debug_LUNDfragmentation
//  G4cout<<"Min mass with Qleft and Qright "<<Qright<<" "<<EstimatedMass<<G4endl;
//  G4cout<<"Number_of_quarks "<<Number_of_quarks<<" Number_of_squarks "<<Number_of_squarks<<G4endl;
#endif

	if(Number_of_quarks==4)
	{
		if(StringM > 1880.) {       // 382. Uzhi May 2016            // 2*Mn = 1880
          	  if(Number_of_squarks==0)      {EstimatedMass += 1320.*MeV;}//560+1320=1880=2*Mn
    	          else if(Number_of_squarks==1) {EstimatedMass += 1150.*MeV;}//920+1150=2070=M(Lam+N)
          	  else if(Number_of_squarks==2) {EstimatedMass +=  960.*MeV;}//1280+960=2240= 2*M Lam
          	  else if(Number_of_squarks==3) {EstimatedMass +=  800.*MeV;}//1640+800=2440=Mxi+Mlam
          	  else if(Number_of_squarks==4) {EstimatedMass +=  640.*MeV;}//2000+640=2640=2*Mxi
                  else {}
                }
/*
		if((StringM > 1880.) && ( EstimatedMass < 2100))     {EstimatedMass = 2020.;}//1880.;}
		//   if((StringM > 1880.) && ( EstimatedMass < 2100))     {EstimatedMass = 2051.;}
		else if((StringM > 2232.) && ( EstimatedMass < 2730)){EstimatedMass = 2570.;}
		else if((StringM > 5130.) && ( EstimatedMass < 3450)){EstimatedMass = 5130.;}
*/
		else
		{
          	  if(Number_of_squarks < 3)     {EstimatedMass -= 200.*MeV;}
          	  else if(Number_of_squarks==3) {EstimatedMass -=  50.*MeV;}
          	  else if(Number_of_squarks==4) {EstimatedMass -=  40.*MeV;}
                  else {}
		}
	}

#ifdef debug_LUNDfragmentation
//  G4cout<<"EstimatedMass -------------------- "<<EstimatedMass <<G4endl;
#endif

	MinimalStringMass=EstimatedMass;
	SetMinimalStringMass2(EstimatedMass);
}

//--------------------------------------------------------------------------------------
void G4LundStringFragmentation::SetMinimalStringMass2(const G4double aValue)
{
	MinimalStringMass2=aValue * aValue;
}
