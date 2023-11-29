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
#include "G4LundStringFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4FragmentingString.hh"
#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4HadronicParameters.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

//#define debug_LUNDfragmentation

// Class G4LundStringFragmentation 
//*************************************************************************************

G4LundStringFragmentation::G4LundStringFragmentation()
  : G4VLongitudinalStringDecay("LundStringFragmentation")
{
    SetMassCut(210.*MeV);   //  Mpi + Delta
                            // For ProduceOneHadron it is required
                            // that no one pi-meson can be produced.
    SigmaQT = 0.435 * GeV;
    Tmt = 190.0 * MeV;       

    SetStringTensionParameter(1.*GeV/fermi);
    SetDiquarkBreakProbability(0.3);

    SetStrangenessSuppression((1.0 - 0.12)/2.0);
    SetDiquarkSuppression(0.07);

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
    
    SetMinMasses();  // For treating of small string decays
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
                      <<"4Mom "<<theString.Get4Momentum()<<G4endl
                      <<"------------------------------------"<<G4endl;
        G4cout<<"String ends Direct "<<theString.GetLeftParton()->GetPDGcode()<<" "
                                     <<theString.GetRightParton()->GetPDGcode()<<" "
                                     <<theString.GetDirection()<< G4endl;
        G4cout<<"Left  mom "<<theString.GetLeftParton()->Get4Momentum()<<G4endl;
        G4cout<<"Right mom "<<theString.GetRightParton()->Get4Momentum()<<G4endl<<G4endl;
        G4cout<<"Check for Fragmentation "<<G4endl;
        #endif

	G4KineticTrackVector * LeftVector(0);

	if (!aString.IsAFourQuarkString() && !IsItFragmentable(&aString))
	{
                #ifdef debug_LUNDfragmentation
                G4cout<<"Non fragmentable - the string is converted to one hadron "<<G4endl;
                #endif
                // SetMassCut(210.*MeV);  // For ProduceOneHadron it is required
                                          // that no one pi-meson can be produced.

		G4double Mcut = GetMassCut();
		SetMassCut(10000.*MeV);
		LeftVector=ProduceOneHadron(&theString);
		SetMassCut(Mcut);

		if ( LeftVector )
		{
		  if ( LeftVector->size() > 0)
                  {
		        LeftVector->operator[](0)->SetFormationTime(theString.GetTimeOfCreation());
		        LeftVector->operator[](0)->SetPosition(theString.GetPosition());
                  }
		  if (LeftVector->size() > 1)
                  {
		        // 2 hadrons created from qq-qqbar are stored
			LeftVector->operator[](1)->SetFormationTime(theString.GetTimeOfCreation());
			LeftVector->operator[](1)->SetPosition(theString.GetPosition());
		  }
		}		
		return LeftVector;
	}

        #ifdef debug_LUNDfragmentation
        G4cout<<"The string will be fragmented. "<<G4endl;
        #endif

	// The string can fragment. At least two particles can be produced.
			       LeftVector =new G4KineticTrackVector;
	G4KineticTrackVector * RightVector=new G4KineticTrackVector;

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
	while (!RightVector->empty())
	{
		LeftVector->push_back(RightVector->back());
		RightVector->erase(RightVector->end()-1);
	}
	delete RightVector;

	return LeftVector;
}

//----------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::IsItFragmentable(const G4FragmentingString * const string)
{
	SetMinimalStringMass(string);
        //G4cout<<"MinM StrM "<<MinimalStringMass<<" "<< string->Get4Momentum().mag()<<G4endl;

	return std::abs(MinimalStringMass) < string->Get4Momentum().mag();

        //MinimalStringMass is negative and large for a string with unknown particles in a final 2-particle decay.
}

//----------------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::Loop_toFragmentString( const G4ExcitedString  &theString,
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

        G4LorentzRotation toCmsI, toObserverFrameI;

	G4bool final_success=false;
	G4bool inner_success=true;

	G4int attempt=0;

	while ( ! final_success && attempt++ < StringLoopInterrupt )
	{       // If the string fragmentation does not be happend, 
	        // repeat the fragmentation.

                G4FragmentingString *currentString = new G4FragmentingString(theString);
                toCmsI = currentString->TransformToAlignedCms();
                toObserverFrameI = toCmsI.inverse();

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
                        G4cout<<"The string will fragment. "<<G4endl;;
                        //G4cout<<"1 "<<currentString->GetDecayDirection()<<G4endl;
                        #endif
			G4FragmentingString *newString=0;  // used as output from SplitUp.

			toCms=currentString->TransformToAlignedCms();
                        toObserverFrame= toCms.inverse();
			
                        #ifdef debug_LUNDfragmentation
                        //G4cout<<"CMS Left  mom "<<currentString->GetPleft()<<G4endl;
                        //G4cout<<"CMS Right mom "<<currentString->GetPright()<<G4endl;
                        //G4cout<<"CMS String M  "<<currentString->GetPstring()<<G4endl;
                        #endif

			G4KineticTrack * Hadron=Splitup(currentString,newString);

			if ( Hadron != 0 )  // Store the hadron                               
			{
                                #ifdef debug_LUNDfragmentation
                                G4cout<<"Hadron prod at fragm. "<<Hadron->GetDefinition()->GetParticleName()<<G4endl;
                                //G4cout<<"2 "<<currentString->GetDecayDirection()<<G4endl;
                                #endif

				Hadron->Set4Momentum(toObserverFrame*Hadron->Get4Momentum());

				G4double TimeOftheStringCreation=theString.GetTimeOfCreation();
				G4ThreeVector PositionOftheStringCreation(theString.GetPosition());

				G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
				G4LorentzVector Momentum = toObserverFrame*Coordinate;
				Hadron->SetFormationTime(TimeOftheStringCreation + Momentum.e() - fermi/c_light);
				G4ThreeVector aPosition(Momentum.vect());
				Hadron->SetPosition(PositionOftheStringCreation+aPosition);
 
                                // Open to protect hadron production at fragmentation 
				if ( currentString->GetDecayDirection() > 0 )
                                {
					LeftVector->push_back(Hadron); 
                                } else
                                {
					RightVector->push_back(Hadron);
                                }
				delete currentString;
				currentString=newString;
			} else {
                          if ( newString ) delete newString;
                        }

                        currentString->LorentzRotate(toObserverFrame);
		};

                if ( loopCounter >= maxNumberOfLoops ) {
                  inner_success=false;
                }

		// Split remaining string into 2 final hadrons.
                #ifdef debug_LUNDfragmentation
                if (inner_success) G4cout<<"Split remaining string into 2 final hadrons."<<G4endl;
                #endif

		if ( inner_success && SplitLast(currentString, LeftVector, RightVector) )  // Close to protect Last Str. Decay
		{
		  final_success = true;
		}

		delete currentString;
	}  // End of the loop where we try to fragment the string.

        G4int sign = +1;
        if ( theString.GetDirection() < 0 ) sign = -1;
        for ( unsigned int hadronI = 0; hadronI < LeftVector->size(); ++hadronI ) {
           G4LorentzVector Tmp = LeftVector->operator[](hadronI)->Get4Momentum();
           Tmp.setZ(sign*Tmp.getZ());
           Tmp *= toObserverFrameI;
           LeftVector->operator[](hadronI)->Set4Momentum(Tmp);
        }
        for ( unsigned int hadronI = 0; hadronI < RightVector->size(); ++hadronI ) {
           G4LorentzVector Tmp = RightVector->operator[](hadronI)->Get4Momentum();
           Tmp.setZ(sign*Tmp.getZ());
           Tmp *= toObserverFrameI;
           RightVector->operator[](hadronI)->Set4Momentum(Tmp);
        }

	return final_success;
}

//----------------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::StopFragmenting(const G4FragmentingString * const string)
{
	SetMinimalStringMass(string); 

	if ( MinimalStringMass < 0.) return true;

	if (string->IsAFourQuarkString())
	{
		return G4UniformRand() < G4Exp(-0.0005*(string->Mass() - MinimalStringMass));
	} else {
           
                if (MinimalStringMass < 0.0 ) return false;  // For a string with di-quark having c or b quarks and s, c, b quarks

		G4bool Result = G4UniformRand() < 
				G4Exp(-0.66e-6*(string->Mass()*string->Mass() - MinimalStringMass*MinimalStringMass));
                // G4bool Result = string->Mass() < MinimalStringMass + 150.*MeV*G4UniformRand();     // a'la LUND

                #ifdef debug_LUNDfragmentation 
                G4cout<<"StopFragmenting MinimalStringMass string->Mass() "<<MinimalStringMass
                      <<" "<<string->Mass()<<G4endl;
                G4cout<<"StopFragmenting - Yes/No "<<Result<<G4endl;
                #endif	
		return Result;
	}
}

//-----------------------------------------------------------------------------

G4KineticTrack * G4LundStringFragmentation::Splitup(G4FragmentingString *string, 
			                            G4FragmentingString *&newString)
{
       #ifdef debug_LUNDfragmentation
       G4cout<<G4endl;
       G4cout<<"Start SplitUP ========================="<<G4endl;
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

       G4double StringMass=string->Mass();

       G4double ProbDqADq = GetDiquarkSuppress();
       G4double ProbSaS   = 1.0 - 2.0 * GetStrangeSuppress();

       #ifdef debug_LUNDfragmentation
       G4cout<<"StrMass DiquarkSuppression           "<<StringMass<<" "<<GetDiquarkSuppress()<<G4endl; 
       #endif   

       G4int NumberOfpossibleBaryons = 2;

       if (string->GetLeftParton()->GetParticleSubType()  != "quark") NumberOfpossibleBaryons++; 
       if (string->GetRightParton()->GetParticleSubType() != "quark") NumberOfpossibleBaryons++;

       G4double ActualProb  = ProbDqADq ;
       ActualProb *= (1.0-G4Pow::GetInstance()->powA(NumberOfpossibleBaryons*1400.0/StringMass, 8.0));
       if(ActualProb <0.0) ActualProb = 0.;

       SetDiquarkSuppression(ActualProb); 

       G4double Mth = 1250.0;                                  // 2 Mk + Mpi
       if ( NumberOfpossibleBaryons == 3 ){Mth = 2520.0;}       // Mlambda/Msigma + Mk + Mpi
       else if ( NumberOfpossibleBaryons == 4 ){Mth = 2380.0;}  // 2 Mlambda/Msigma + Mk + Mpi
       else {}

       ActualProb = ProbSaS;
       ActualProb *= (1.0 - G4Pow::GetInstance()->powA( Mth/StringMass, 2.5 ));
       if ( ActualProb < 0.0 ) ActualProb = 0.0;
       SetStrangenessSuppression((1.0-ActualProb)/2.0);

       #ifdef debug_LUNDfragmentation
       G4cout<<"StrMass DiquarkSuppression corrected "<<StringMass<<" "<<GetDiquarkSuppress()<<G4endl; 
       #endif

       if (string->DecayIsQuark())
       {
          HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);
       } else {
          HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
       }

       SetDiquarkSuppression(ProbDqADq);
       SetStrangenessSuppression((1.0-ProbSaS)/2.0);

       if ( HadronDefinition == NULL ) { G4KineticTrack * Hadron =0; return Hadron; }

       #ifdef debug_LUNDfragmentation
       G4cout<<"The parton "<<string->GetDecayParton()->GetPDGEncoding()<<" "
             <<" produces hadron "<<HadronDefinition->GetParticleName()
             <<" and is transformed to "<<newStringEnd->GetPDGEncoding()<<G4endl;
       G4cout<<"The side of the string decay Left/Right (1/-1) "<<SideOfDecay<<G4endl;
       #endif
       // create new String from old, ie. keep Left and Right order, but replace decay

       if ( newString ) delete newString;

       newString=new G4FragmentingString(*string,newStringEnd);  // To store possible quark containt of new string

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
	   
           if ( newString ) delete newString;

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

//-----------------------------------------------------------------------------

G4ParticleDefinition * G4LundStringFragmentation::DiQuarkSplitup(G4ParticleDefinition* decay,
                                                                 G4ParticleDefinition *&created)
{
   G4double StrSup=GetStrangeSuppress();
   G4double ProbQQbar = (1.0 - 2.0*StrSup)*1.25;

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

      G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1;  // if we have a quark, we need antiquark

      SetStrangenessSuppression((1.0-ProbQQbar)/2.0);
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      SetStrangenessSuppression((1.0-StrSup)/2.0);

      //... Build new Diquark
      G4int QuarkEncoding=QuarkPair.second->GetPDGEncoding();
      G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
      G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
      created = FindParticle(NewDecayEncoding);
      G4ParticleDefinition * decayQuark=FindParticle(decayQuarkEncoding);
      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decayQuark);
      StrangeSuppress=StrSup;

      return had;

   } else {
      //... Diquark does not break

      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1;  // if we have a diquark, we need quark

      StrangeSuppress=(1.0 - ProbQQbar)/2.0;
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted

      created = QuarkPair.second;

      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
      StrangeSuppress=StrSup;

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

       if ( MinimalStringMass < 0.0 ) return nullptr;

        #ifdef debug_LUNDfragmentation
        G4cout<<G4endl<<"Start LUND SplitEandP "<<G4endl;
        G4cout<<"String 4 mom, String M and Mt "<<String4Momentum<<" "<<String4Momentum.mag()
              <<" "<<std::sqrt(StringMT2)<<G4endl;
        G4cout<<"Hadron "<<pHadron->GetParticleName()<<G4endl;
        G4cout<<"HadM MinimalStringMassLeft StringM hM+sM "<<HadronMass<<" "<<MinimalStringMass<<" "
              <<String4Momentum.mag()<<" "<<HadronMass+MinimalStringMass<<G4endl;
        #endif

	if ((HadronMass + MinimalStringMass > string->Mass()) || MinimalStringMass < 0.) 
	{
          #ifdef debug_LUNDfragmentation
          G4cout<<"Mass of the string is not sufficient to produce the hadron!"<<G4endl;
          #endif
	  return 0;
	}  // have to start all over!

	String4Momentum.setPz(0.);
	G4ThreeVector StringPt=String4Momentum.vect();
        StringPt.setZ(0.);
	
	// calculate and assign hadron transverse momentum component HadronPx and HadronPy
	G4ThreeVector HadronPt    , RemSysPt; 
	G4double      HadronMassT2, ResidualMassT2;
	G4double HadronMt, Pt, Pt2, phi;

        G4double TmtCur = Tmt;

        if ( (string->GetDecayParton()->GetParticleSubType()== "quark") &&
             (pHadron->GetBaryonNumber() != 0) ) {
          TmtCur = Tmt*0.37;			        // q->B     
        } else if ( (string->GetDecayParton()->GetParticleSubType()== "quark") &&
                    (pHadron->GetBaryonNumber() == 0) ) {
          //TmtCur = Tmt;                               // q->M
	} else if ( (string->GetDecayParton()->GetParticleSubType()== "di_quark") && 
                    (pHadron->GetBaryonNumber() == 0) ) {
          //TmtCur = Tmt*0.89;                          // qq -> M
        } else if ( (string->GetDecayParton()->GetParticleSubType()== "di_quark") &&
                    (pHadron->GetBaryonNumber() != 0) ) {
          TmtCur = Tmt*1.35;                            // qq -> B
        }

        //...  sample Pt of the hadron
        G4int attempt=0;
        do
        {
          attempt++; if (attempt > StringLoopInterrupt) {return 0;}

          HadronMt = HadronMass - TmtCur*G4Log(G4UniformRand());
	  Pt2 = sqr(HadronMt)-sqr(HadronMass); Pt=std::sqrt(Pt2);
	  phi = 2.*pi*G4UniformRand();
          HadronPt = G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0. );
          RemSysPt = StringPt - HadronPt;
          HadronMassT2 = sqr(HadronMass) + HadronPt.mag2();
          ResidualMassT2=sqr(MinimalStringMass) + RemSysPt.mag2();

        } while (std::sqrt(HadronMassT2) + std::sqrt(ResidualMassT2) > StringMT);

	//...  sample z to define hadron longitudinal momentum and energy
	//... but first check the available phase space

	G4double Pz2 = (sqr(StringMT2 - HadronMassT2 - ResidualMassT2) -
			4*HadronMassT2 * ResidualMassT2)/4./StringMT2;

	if (Pz2 < 0 ) {return 0;}          // have to start all over!

	//... then compute allowed z region  z_min <= z <= z_max

	G4double Pz = std::sqrt(Pz2);
	G4double zMin = (std::sqrt(HadronMassT2+Pz2) - Pz)/std::sqrt(StringMT2);
        // G4double zMin = (std::sqrt(HadronMassT2+Pz2) - 0.)/std::sqrt(StringMT2); // For testing purposes
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

	if (!((std::abs(PDGEncodingOfDecayParton) > 1000) && (HadronEncoding > 1000)) )
	{    // ---------------- Quark fragmentation  and qq-> meson ----------------------
           Alund=1.;
           Blund=0.7/GeV/GeV;

	   G4double BMt2 = Blund*Mt2;
	   if (Alund == 1.0) {
	     zOfMaxyf=BMt2/(Blund*Mt2 + 1.);}
	   else {  
   	     zOfMaxyf = ((1.0+BMt2) - std::sqrt(sqr(1.0-BMt2) + 4.0*BMt2*Alund))/2.0/(1.-Alund);
	   }

	   if (zOfMaxyf < zmin) {zOfMaxyf=zmin;}
	   if (zOfMaxyf > zmax) {zOfMaxyf=zmax;}
	   maxYf=(1-zOfMaxyf)/zOfMaxyf * G4Exp(-Blund*Mt2/zOfMaxyf);

           const G4int maxNumberOfLoops = 1000;
           G4int loopCounter = 0;
	   do
	   {
		z = zmin + G4UniformRand()*(zmax-zmin);
                //yf = (1-z)/z * G4Exp(-Blund*Mt2/z);
		yf = G4Pow::GetInstance()->powA(1.0-z,Alund)/z*G4Exp(-BMt2/z);
	   }
	   while ( (G4UniformRand()*maxYf > yf) && ++loopCounter < maxNumberOfLoops );
           if ( loopCounter >= maxNumberOfLoops ) {
             z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
           }
	   return z;
        }

	if (std::abs(PDGEncodingOfDecayParton) > 1000)     
	{
                G4double an = 2.5;
                an +=(sqr(Px)+sqr(Py))/sqr(GeV)-0.5;
                z=zmin + (zmax-zmin)*G4Pow::GetInstance()->powA(G4UniformRand(),1./an);
                if( PDGEncodingOfDecayParton > 3000 ) z=zmin+zmax-z;
	}

	return z;
}

//----------------------------------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::SplitLast(G4FragmentingString * string,
                                            G4KineticTrackVector * LeftVector,
                                            G4KineticTrackVector * RightVector)
{
	//... perform last cluster decay
        SetMinimalStringMass( string);
        if ( MinimalStringMass < 0.) return false;
        #ifdef debug_LUNDfragmentation
        G4cout<<G4endl<<"Split last-----------------------------------------"<<G4endl;
        G4cout<<"MinimalStringMass "<<MinimalStringMass<<G4endl;
        G4cout<<"Left  "<<string->GetLeftParton()->GetPDGEncoding()<<" "<<string->GetPleft()<<G4endl;
        G4cout<<"Right "<<string->GetRightParton()->GetPDGEncoding()<<" "<<string->GetPright()<<G4endl;
        G4cout<<"String4mom "<<string->GetPstring()<<" "<<string->GetPstring().mag()<<G4endl;
        #endif

        G4LorentzVector Str4Mom=string->Get4Momentum();
        G4LorentzRotation toCms(-1*Str4Mom.boostVector());
        G4LorentzVector Pleft = toCms * string->GetPleft();
        toCms.rotateZ(-1*Pleft.phi());
        toCms.rotateY(-1*Pleft.theta());
	
        G4LorentzRotation toObserverFrame= toCms.inverse();

	G4double StringMass=string->Mass();

	G4ParticleDefinition * LeftHadron(0), * RightHadron(0);

        NumberOf_FS=0;
	for (G4int i=0; i<350; i++) {FS_Weight[i]=0.;}

	G4int sampledState = 0;

        #ifdef debug_LUNDfragmentation
        G4cout<<"StrMass "<<StringMass<<" q's "
              <<string->GetLeftParton()->GetParticleName()<<" "
              <<string->GetRightParton()->GetParticleName()<<G4endl;
        #endif

	string->SetLeftPartonStable(); // to query quark contents..

	if (string->IsAFourQuarkString() )
	{
          G4int IDleft =std::abs(string->GetLeftParton()->GetPDGEncoding());
          G4int IDright=std::abs(string->GetRightParton()->GetPDGEncoding());

          if ( (IDleft > 3000) || (IDright > 3000) ) {
            if ( ! Diquark_AntiDiquark_belowThreshold_lastSplitting(string, LeftHadron, RightHadron) )
            {
              return false;
            } 
          } else {
		// The string is qq-qqbar type. Diquarks are on the string ends
	        if (StringMass-MinimalStringMass < 0.)
		{
			if (! Diquark_AntiDiquark_belowThreshold_lastSplitting(string, LeftHadron, RightHadron) )
                        {
				return false;
                        }
		} else
		{
			Diquark_AntiDiquark_aboveThreshold_lastSplitting(string, LeftHadron, RightHadron);

			if (NumberOf_FS == 0) return false;

                        sampledState = SampleState();
			if (string->GetLeftParton()->GetPDGEncoding() < 0)
			{
				LeftHadron =FS_LeftHadron[sampledState];
				RightHadron=FS_RightHadron[sampledState];
			} else
			{
				LeftHadron =FS_RightHadron[sampledState];
				RightHadron=FS_LeftHadron[sampledState];
			}
		}
          }  // ID > 3300
        } else {
		if (string->DecayIsQuark() && string->StableIsQuark() )
		{       //... there are quarks on cluster ends
                        #ifdef debug_LUNDfragmentation
                        G4cout<<"Q Q string LastSplit"<<G4endl;
                        #endif

			Quark_AntiQuark_lastSplitting(string, LeftHadron, RightHadron);
		
			if (NumberOf_FS == 0) return false;
                	sampledState = SampleState();
			if (string->GetLeftParton()->GetPDGEncoding() < 0)
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

			if (NumberOf_FS == 0) return false;                           
                	sampledState = SampleState();

			if (string->GetLeftParton()->GetParticleSubType() == "quark")
			{
				LeftHadron =FS_LeftHadron[sampledState];
				RightHadron=FS_RightHadron[sampledState];
			} else 
			{
				LeftHadron =FS_RightHadron[sampledState];
				RightHadron=FS_LeftHadron[sampledState];
			}			
		}

	}

        #ifdef debug_LUNDfragmentation
        G4cout<<"Sampled hadrons: "<<LeftHadron->GetParticleName()<<" "<<RightHadron->GetParticleName()<<G4endl;
        #endif

	G4LorentzVector P_left  =string->GetPleft(), P_right = string->GetPright();

	G4LorentzVector  LeftMom, RightMom;
	G4ThreeVector    Pos;

	Sample4Momentum(&LeftMom,  LeftHadron->GetPDGMass(),
			&RightMom, RightHadron->GetPDGMass(),
			StringMass);

        // Sample4Momentum ascribes LeftMom.pz() along positive Z axis for baryons in many cases.
        // It must be negative in case when the system is moving against Z axis.

	if (!(string->DecayIsQuark() && string->StableIsQuark() ))
	{ // Only for qq - q, q - qq, and qq - qqbar -------------------

          if ( G4UniformRand() <= 0.5 )
	  {
	    if (P_left.z() <= 0.) {G4LorentzVector tmp = LeftMom; LeftMom=RightMom; RightMom=tmp;}
	  } 
	  else  
	  {
	    if (P_right.z() >= 0.) {G4LorentzVector tmp = LeftMom; LeftMom=RightMom; RightMom=tmp;}
	  }
	}

	LeftMom *=toObserverFrame;
	RightMom*=toObserverFrame;

	LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
	RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

	string->LorentzRotate(toObserverFrame);
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
        G4bool isOK = false;
	do
	{
		G4int LeftQuark1= string->GetLeftParton()->GetPDGEncoding()/1000;
		G4int LeftQuark2=(string->GetLeftParton()->GetPDGEncoding()/100)%10;

		G4int RightQuark1= string->GetRightParton()->GetPDGEncoding()/1000;
		G4int RightQuark2=(string->GetRightParton()->GetPDGEncoding()/100)%10;

		if (G4UniformRand()<0.5)
		{
			LeftHadron =hadronizer->Build(FindParticle( LeftQuark1),
						      FindParticle(RightQuark1));
			RightHadron= (LeftHadron == nullptr) ? nullptr :
                                                      hadronizer->Build(FindParticle( LeftQuark2),
						      FindParticle(RightQuark2));
		} else
		{
			LeftHadron =hadronizer->Build(FindParticle( LeftQuark1),
						      FindParticle(RightQuark2));
			RightHadron=(LeftHadron == nullptr) ? nullptr :
			                              hadronizer->Build(FindParticle( LeftQuark2),
						      FindParticle(RightQuark1));
		}

		isOK = (LeftHadron != nullptr) && (RightHadron != nullptr);

		if(isOK) { isOK = (StringMass > LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()); }
		++cClusterInterrupt;
		//... repeat procedure, if mass of cluster is too low to produce hadrons
		//... ClusterMassCut = 0.15*GeV model parameter
	}
	while (isOK == false && cClusterInterrupt < ClusterLoopInterrupt);
	/* Loop checking, 07.08.2015, A.Ribon */
  	return isOK;
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

	if (string->GetLeftParton()->GetPDGEncoding() < 0)
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
	for (G4int ProdQ=1; ProdQ < 6; ProdQ++)
	{
		G4int StateADiQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do  // while(Meson[AbsIDquark-1][ProdQ-1][StateQ]<>0);
		{
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(
							-Baryon[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]);

			if (LeftHadron == NULL) continue;
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			G4int StateDiQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do // while(Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]<>0);
			{
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(
								+Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);

				if (RightHadron == NULL) continue;
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if (StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 349 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Diquark_AntiDiquark_aboveThreshold_lastSplitting ",
                                                       "HAD_LUND_001", JustWarning, ed );
                                          NumberOf_FS = 349;
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
				} // End of if (StringMass > LeftHadronMass + RightHadronMass)

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
	} // End of for (G4int ProdQ=1; ProdQ < 4; ProdQ++)

  	return true;
}

//----------------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::Quark_Diquark_lastSplitting(G4FragmentingString * & string,
                                                              G4ParticleDefinition * & LeftHadron,
                                                              G4ParticleDefinition * & RightHadron)
{
	G4double StringMass   = string->Mass();
	G4double StringMassSqr= sqr(StringMass);

	G4ParticleDefinition * Di_Quark;
	G4ParticleDefinition * Quark;

	if (string->GetLeftParton()->GetParticleSubType()== "quark")
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
	G4int SignDiQ= 1;
	if (IDdi_quark < 0) SignDiQ=-1;

	NumberOf_FS=0;
	for (G4int ProdQ=1; ProdQ < 4; ProdQ++)  // Loop over quark-antiquark cases: u-ubar, d-dbar, s-sbar 
	{                                        // (as last splitting, do not consider c-cbar and b-bbar cases)
		G4int SignQ;
		if (IDquark > 0)
		{
		        SignQ=-1;
		        if (IDquark == 2)                   SignQ= 1;
			if ((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0
			if ((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar
		        if (IDquark == 4)                   SignQ= 1; // D+, D0, Ds+
			if (IDquark == 5)                   SignQ=-1; // B-, anti_B0, anti_Bs0
		} else
		{
			SignQ= 1;
			if (IDquark == -2)                  SignQ=-1;
			if ((IDquark ==-1) && (ProdQ == 3)) SignQ=-1; // K0bar
			if ((IDquark ==-3) && (ProdQ == 1)) SignQ= 1; // K0
			if (IDquark == -4)                  SignQ=-1; // D-, anti_D0, anti_Ds+
			if (IDquark == -5)                  SignQ= 1; // B+, B0, Bs0
		}

		if (AbsIDquark == ProdQ)            SignQ= 1;

		G4int StateQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do  // while(Meson[AbsIDquark-1][ProdQ-1][StateQ]<>0);
		{
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignQ*
							Meson[AbsIDquark-1][ProdQ-1][StateQ]);
			if (LeftHadron == NULL) continue;
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			G4int StateDiQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do // while(Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]<>0);
			{
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignDiQ*
								Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);
				if (RightHadron == NULL) continue;
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if (StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 349 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Quark_Diquark_lastSplitting ",
                                                       "HAD_LUND_002", JustWarning, ed );
                                          NumberOf_FS = 349;
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
				} // End of if (StringMass > LeftHadronMass + RightHadronMass)

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
	}

	return true;
}

//----------------------------------------------------------------------------------------

G4bool G4LundStringFragmentation::Quark_AntiQuark_lastSplitting(G4FragmentingString * & string,
                                                                G4ParticleDefinition * & LeftHadron,
                                                                G4ParticleDefinition * & RightHadron)
{
	G4double StringMass   = string->Mass();
	G4double StringMassSqr= sqr(StringMass);

	G4ParticleDefinition * Quark;
	G4ParticleDefinition * Anti_Quark;

	if (string->GetLeftParton()->GetPDGEncoding()>0)
	{
		Quark     =string->GetLeftParton();
		Anti_Quark=string->GetRightParton();
	} else
	{
		Quark     =string->GetRightParton();
		Anti_Quark=string->GetLeftParton();
	}

	G4int IDquark         =Quark->GetPDGEncoding();
	G4int AbsIDquark      =std::abs(IDquark);
        G4int QuarkCharge     =Qcharge[IDquark-1];

	G4int IDanti_quark    =Anti_Quark->GetPDGEncoding();
	G4int AbsIDanti_quark =std::abs(IDanti_quark);
        G4int AntiQuarkCharge =-Qcharge[AbsIDanti_quark-1];

        G4int LeftHadronCharge(0), RightHadronCharge(0);

        //G4cout<<"Q Qbar "<<IDquark<<" "<<IDanti_quark<<G4endl;

	NumberOf_FS=0;
	for (G4int ProdQ=1; ProdQ < 4; ProdQ++)  // Loop over quark-antiquark cases: u-ubar, d-dbar, s-sbar 
	{                                        // (as last splitting, do not consider c-cbar and b-bbar cases)
                //G4cout<<"NumberOf_FS ProdQ "<<NumberOf_FS<<" "<<ProdQ<<G4endl;
		LeftHadronCharge = QuarkCharge - Qcharge[ProdQ-1];
		G4int SignQ = LeftHadronCharge/3; if (SignQ == 0) SignQ = 1;

		if ((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0       (d,sbar)
		if ((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar    (s,dbar)
                if ((IDquark == 4) && (ProdQ == 2)) SignQ= 1; // D0       (c,ubar)
                if ((IDquark == 5) && (ProdQ == 1)) SignQ=-1; // anti_B0  (b,dbar)
                if ((IDquark == 5) && (ProdQ == 3)) SignQ=-1; // anti_Bs0 (b,sbar)
		
                RightHadronCharge = AntiQuarkCharge + Qcharge[ProdQ-1];
		G4int SignAQ = RightHadronCharge/3; if (SignAQ == 0) SignAQ = 1;

		if ((IDanti_quark ==-1) && (ProdQ == 3)) SignAQ=-1; // K0bar   (dbar,s)
		if ((IDanti_quark ==-3) && (ProdQ == 1)) SignAQ= 1; // K0      (sbar,d)
		if ((IDanti_quark ==-4) && (ProdQ == 2)) SignAQ=-1; // anti_D0 (cbar,u)
		if ((IDanti_quark ==-5) && (ProdQ == 1)) SignAQ= 1; // B0      (bbar,d)
		if ((IDanti_quark ==-5) && (ProdQ == 3)) SignAQ= 1; // Bs0     (bbar,s)

                //G4cout<<"ProQ signs "<<ProdQ<<" "<<SignQ<<" "<<SignAQ<<G4endl;

		G4int StateQ=0;
                const G4int maxNumberOfLoops = 1000;
                G4int loopCounter = 0;
		do
		{
                        //G4cout<<"[AbsIDquark-1][ProdQ-1][StateQ "<<AbsIDquark-1<<" "
                        //<<ProdQ-1<<" "<<StateQ<<" "<<SignQ*Meson[AbsIDquark-1][ProdQ-1][StateQ]<<G4endl;
			LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignQ*
						       Meson[AbsIDquark-1][ProdQ-1][StateQ]);
                        //G4cout<<"LeftHadron "<<LeftHadron<<G4endl; 
			if (LeftHadron == NULL) { StateQ++; continue; }
                        //G4cout<<"LeftHadron "<<LeftHadron->GetParticleName()<<G4endl;
			G4double LeftHadronMass=LeftHadron->GetPDGMass();

			G4int StateAQ=0;
                        const G4int maxNumberOfInternalLoops = 1000;
                        G4int internalLoopCounter = 0;
			do
			{
                                //G4cout<<"           [AbsIDanti_quark-1][ProdQ-1][StateAQ] "<<AbsIDanti_quark-1<<" "
                                //      <<ProdQ-1<<" "<<StateAQ<<" "<<SignAQ*Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]<<G4endl;
				RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(SignAQ*
								Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]);
                                //G4cout<<"RightHadron "<<RightHadron<<G4endl;
				if(RightHadron == NULL) { StateAQ++; continue; }
                                //G4cout<<"RightHadron "<<RightHadron->GetParticleName()<<G4endl;
				G4double RightHadronMass=RightHadron->GetPDGMass();

				if (StringMass > LeftHadronMass + RightHadronMass)
				{
                                        if ( NumberOf_FS > 349 ) {
                                          G4ExceptionDescription ed;
                                          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
                                          G4Exception( "G4LundStringFragmentation::Quark_AntiQuark_lastSplitting ",
                                                       "HAD_LUND_003", JustWarning, ed );
                                          NumberOf_FS = 349;
                                        }

                                        G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),
								sqr(RightHadronMass));
					//FS_Psqr=1.;
					FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*
							       MesonWeight[AbsIDquark-1][ProdQ-1][StateQ]*
							       MesonWeight[AbsIDanti_quark-1][ProdQ-1][StateAQ]*
							       Prob_QQbar[ProdQ-1];
					if (string->GetLeftParton()->GetPDGEncoding()>0)
					{
						FS_LeftHadron[NumberOf_FS] = RightHadron;
						FS_RightHadron[NumberOf_FS]= LeftHadron;
					} else
					{
						FS_LeftHadron[NumberOf_FS] = LeftHadron;
						FS_RightHadron[NumberOf_FS]= RightHadron;
					}

					NumberOf_FS++;

				}

				StateAQ++;
                                //G4cout<<"               StateAQ Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ] "<<StateAQ<<" "
                                //      <<Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]<<" "<<internalLoopCounter<<G4endl;
			} while ( (Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]!=0) &&
                                  ++internalLoopCounter < maxNumberOfInternalLoops );
                          if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
                            return false;
                          }

			StateQ++;
                        //G4cout<<"StateQ Meson[AbsIDquark-1][ProdQ-1][StateQ] "<<StateQ<<" "
                        //      <<Meson[AbsIDquark-1][ProdQ-1][StateQ]<<" "<<loopCounter<<G4endl;

		} while ( (Meson[AbsIDquark-1][ProdQ-1][StateQ]!=0) &&
                         ++loopCounter < maxNumberOfLoops );
                  if ( loopCounter >= maxNumberOfLoops ) {
                    return false;
                  }
	} // End of for (G4int ProdQ=1; ProdQ < 4; ProdQ++)

	return true;
}

//----------------------------------------------------------------------------------------------------------

G4int G4LundStringFragmentation::SampleState(void) 
{
        if ( NumberOf_FS > 349 ) {
          G4ExceptionDescription ed;
          ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
          G4Exception( "G4LundStringFragmentation::SampleState ", "HAD_LUND_004", JustWarning, ed );
          NumberOf_FS = 349;
        }

	G4double SumWeights=0.;
	for (G4int i=0; i<NumberOf_FS; i++) {SumWeights+=FS_Weight[i];}

	G4double ksi=G4UniformRand();
	G4double Sum=0.;
	G4int indexPosition = 0;

	for (G4int i=0; i<NumberOf_FS; i++)
	{
		Sum+=(FS_Weight[i]/SumWeights);
		indexPosition=i;
		if (Sum >= ksi) break;
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
        #ifdef debug_LUNDfragmentation
        G4cout<<"Sampling of momenta of 2 last produced hadrons ----------------"<<G4endl;
        G4cout<<"Init Mass "<<InitialMass<<" FirstM "<<Mass<<" SecondM "<<AntiMass<<" ProbIsotropy "<<G4endl;  
        #endif

	G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) -   
 			sqr(2.*Mass*AntiMass);
	G4double Pabs = (r_val > 0.)? std::sqrt(r_val)/(2.*InitialMass) : 0;

        const G4int maxNumberOfLoops = 1000;
	G4double SigmaQTw = SigmaQT;
	if ( Mass > 930. || AntiMass > 930. ) {
	  SigmaQT *= ( 1.0 - 0.55*sqr( (Mass+AntiMass)/InitialMass ) );
        }
        if ( Mass < 930. && AntiMass < 930. ) {}     // q-qbar string
        if ( ( Mass < 930. && AntiMass > 930. ) ||
	     ( Mass > 930. && AntiMass < 930. ) ) {  // q-di_q string
          //SigmaQT = -1.;  // isotropical decay
        }
        if ( Mass > 930. && AntiMass > 930. ) {      // qq-qqbar string
          SigmaQT *= ( 1.0 - 0.55*sqr( (Mass+AntiMass)/InitialMass ) );
        }

        G4int loopCounter = 0;
	do
	{
		Pt=SampleQuarkPt(Pabs); Pt.setZ(0); G4double Pt2=Pt.mag2();
		MassMt    = std::sqrt(    Mass *     Mass + Pt2);
		AntiMassMt= std::sqrt(AntiMass * AntiMass + Pt2);
	}
	while ( (InitialMass < MassMt + AntiMassMt) && ++loopCounter < maxNumberOfLoops ); 

        SigmaQT = SigmaQTw;

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

        #ifdef debug_LUNDfragmentation
        G4cout<<"Fmass Mom "<<Mom->getX()<<" "<<Mom->getY()<<" "<<Mom->getZ()<<" "<<Mom->getT()<<G4endl;
        G4cout<<"Smass Mom "<<AntiMom->getX()<<" "<<AntiMom->getY()<<" "<<AntiMom->getZ()
              <<" "<<AntiMom->getT()<<G4endl;
        #endif
}

//------------------------------------------------------------------------

G4double G4LundStringFragmentation::lambda(G4double S, G4double m1_Sqr, G4double m2_Sqr)
{ 
	G4double lam = sqr(S - m1_Sqr - m2_Sqr) - 4.*m1_Sqr*m2_Sqr;
	return lam;
}

// --------------------------------------------------------------
G4LundStringFragmentation::~G4LundStringFragmentation()
{}

