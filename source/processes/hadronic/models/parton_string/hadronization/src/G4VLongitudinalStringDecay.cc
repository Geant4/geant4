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
// $Id: G4VLongitudinalStringDecay.cc 102717 2017-02-20 10:37:13Z gcosmo $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
//               redesign  Gunter Folger, August/September 2001
// -----------------------------------------------------------------------------
#include "G4VLongitudinalStringDecay.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4FragmentingString.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleChange.hh"
#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4Gluons.hh"

#include "G4Exp.hh"
#include "G4Log.hh"

//------------------------debug switches
//#define debug_VStringDecay

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
   DiquarkSuppress  = 0.07;
   DiquarkBreakProb = 0.1;
   
   //... pspin_meson is probability to create pseudo-scalar meson 
   pspin_meson = 0.5;

   //... pspin_barion is probability to create 1/2 barion 
   pspin_barion = 0.5;

   //... vectorMesonMix[] is quark mixing parameters for vector mesons (Variable spin = 3)
   vectorMesonMix.resize(6);
   vectorMesonMix[0] = 0.5;  // Uzhi May 2016 //AR-20Oct2014 : it was 0.5
   vectorMesonMix[1] = 0.0;
   vectorMesonMix[2] = 0.5;  // Uzhi May 2016 //AR-20Oct2014 : it was 0.5
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
   Kappa = 1.0 * GeV/fermi;


}


G4VLongitudinalStringDecay::~G4VLongitudinalStringDecay()
   {
   delete hadronizer;
   }

//=============================================================================

// Operators

//-----------------------------------------------------------------------------

int G4VLongitudinalStringDecay::operator==(const G4VLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::operator== forbidden");
	return false;
    }

//-------------------------------------------------------------------------------------

int G4VLongitudinalStringDecay::operator!=(const G4VLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::operator!= forbidden");
	return true;
    }

//***********************************************************************************

// For changing Mass Cut used for selection of very small mass strings
void     G4VLongitudinalStringDecay::SetMassCut(G4double aValue){MassCut=aValue;}
G4double G4VLongitudinalStringDecay::GetMassCut(){return MassCut;}
//-----------------------------------------------------------------------------

// For handling a string with very low mass

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

// The string mass is very low ---------------------------
	
	result=new G4KineticTrackVector;
        
	if ( hadrons.second ==0 )
	{
// Substitute string by light hadron, Note that Energy is not conserved here!

#ifdef debug_VStringDecay
	G4cout << "VlongSF Warning replacing string by single hadron (G4VLongitudinalStringDecay)" <<G4endl;
	G4cout << hadrons.first->GetParticleName()<<G4endl
	       << "string .. " << string->Get4Momentum() << " " 
	       << string->Get4Momentum().m() << G4endl;
#endif		

	       G4ThreeVector   Mom3 = string->Get4Momentum().vect();
	       G4LorentzVector Mom(Mom3, 
	       			   std::sqrt(Mom3.mag2() + 
                                             sqr(hadrons.first->GetPDGMass())));
               result->push_back(new G4KineticTrack(hadrons.first, 0, 
                                                  string->GetPosition(),
                                                          Mom));
	} else 
	{
//... string was qq--qqbar type: Build two stable hadrons,

#ifdef debug_VStringDecay
	G4cout << "VlongSF Warning replacing qq-qqbar string by TWO hadrons (G4VLongitudinalStringDecay)" 
	       << hadrons.first->GetParticleName() << " / " 
	       << hadrons.second->GetParticleName()
	       << "string .. " << string->Get4Momentum() << " " 
	       << string->Get4Momentum().m() << G4endl;
#endif		      

	       G4LorentzVector  Mom1, Mom2;
	       Sample4Momentum(&Mom1, hadrons.first->GetPDGMass(), 
			       &Mom2,hadrons.second->GetPDGMass(),
			       string->Get4Momentum().mag());

	       result->push_back(new G4KineticTrack(hadrons.first, 0, 
                                                    string->GetPosition(), 
                                                            Mom1));
	       result->push_back(new G4KineticTrack(hadrons.second, 0, 
                                                    string->GetPosition(), 
                                                    Mom2));

               G4ThreeVector Velocity = string->Get4Momentum().boostVector();
               result->Boost(Velocity);          
	}

	return result;
	
}

//----------------------------------------------------------------------------------------

G4double G4VLongitudinalStringDecay::FragmentationMass(
            const G4FragmentingString *	const string,
		Pcreate build, pDefPair * pdefs       )
{
        G4double mass;
        static G4ThreadLocal G4bool NeedInit(true);
	static G4ThreadLocal std::vector<double> *nomix_G4MT_TLS_ = 0 ; if (!nomix_G4MT_TLS_) nomix_G4MT_TLS_ = new  std::vector<double>  ;  std::vector<double> &nomix = *nomix_G4MT_TLS_;
	static G4ThreadLocal G4HadronBuilder * minMassHadronizer;
	if ( NeedInit ) 
	{
	   NeedInit = false;
	   nomix.resize(6);
	   for ( G4int i=0; i<6 ; i++ ) nomix[i]=0;

//	   minMassHadronizer=new G4HadronBuilder(pspin_meson,pspin_barion,nomix,nomix);
	   minMassHadronizer=hadronizer;
	}

	if ( build==0 ) build=&G4HadronBuilder::BuildLowSpin;

        G4ParticleDefinition *Hadron1, *Hadron2=0;

        if (!string->FourQuarkString() )
        {
           // spin 0 meson or spin 1/2 barion will be built

           Hadron1 = (minMassHadronizer->*build)(string->GetLeftParton(),
			                         string->GetRightParton());

#ifdef debug_VStringDecay
  G4cout<<"Quarks at the string ends "<<string->GetLeftParton()->GetParticleName()<<" "<<string->GetRightParton()->GetParticleName()<<G4endl;
  G4cout<<"(G4VLongitudinalStringDecay) Hadron "<<Hadron1->GetParticleName()<<" "<<Hadron1->GetPDGMass()<<G4endl;
#endif
           mass= (Hadron1)->GetPDGMass();
        } else
        {
           //... string is qq--qqbar: Build two stable hadrons,
           //...    with extra uubar or ddbar quark pair
	   G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
	   if (string->GetLeftParton()->GetPDGEncoding() < 0) iflc = -iflc;

	   //... theSpin = 4; spin 3/2 baryons will be built
	   Hadron1 = (minMassHadronizer->*build)(string->GetLeftParton(),
                                                 FindParticle(iflc)       );
	   Hadron2 = (minMassHadronizer->*build)(string->GetRightParton(),
                                                 FindParticle(-iflc)      );
           mass = (Hadron1)->GetPDGMass() + (Hadron2)->GetPDGMass();
        }
	
	if ( pdefs != 0 ) 
	{ // need to return hadrons as well....
	   pdefs->first  = Hadron1;
	   pdefs->second = Hadron2;
	}
	   
        return mass;
}

//----------------------------------------------------------------------------

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

//*********************************************************************************
//   For decision on continue or stop string fragmentation
//   virtual G4bool StopFragmenting(const G4FragmentingString  * const string)=0;
//   virtual G4bool IsFragmentable(const G4FragmentingString * const string)=0;

//   If a string can not fragment, make last break into 2 hadrons
//   virtual G4bool SplitLast(G4FragmentingString * string, 
//                            G4KineticTrackVector * LeftVector,
//                            G4KineticTrackVector * RightVector)=0;
//-----------------------------------------------------------------------------
//
//   If a string can fragment, do the following
//
//   For transver of a string to its CMS frame
//-----------------------------------------------------------------------------

G4ExcitedString *G4VLongitudinalStringDecay::CopyExcited(const G4ExcitedString & in)
{
	G4Parton *Left=new G4Parton(*in.GetLeftParton());
	G4Parton *Right=new G4Parton(*in.GetRightParton());
	return new G4ExcitedString(Left,Right,in.GetDirection());
}

//-----------------------------------------------------------------------------
/*                                                        // Uzhi 28 June 2016, shifted to FTF and QGSM
G4KineticTrack * G4VLongitudinalStringDecay::Splitup(
		        G4FragmentingString *string, 
			G4FragmentingString *&newString)
{
G4cout<<"G4VLongStrDeca Splitup ----------------------"<<G4endl;
#ifdef debug_VStringDecay
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
       	   HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);
       } else {
           HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
       }      

#ifdef debug_VStringDecay
  G4cout<<"The parton "<<string->GetDecayParton()->GetPDGEncoding()<<" "
        <<" produces hadron "<<HadronDefinition->GetParticleName()
        <<" and is transformed to "<<newStringEnd->GetPDGEncoding()<<G4endl;
  G4cout<<"The side of the string decay Left/Right (1/-1) "<<SideOfDecay<<G4endl;
#endif
// create new String from old, ie. keep Left and Right order, but replace decay

       newString=new G4FragmentingString(*string,newStringEnd); // To store possible
                                                                // quark containt of new string

#ifdef debug_VStringDecay
  G4cout<<"An attempt to determine its energy (SplitEandP)"<<G4endl;
#endif
       G4LorentzVector* HadronMomentum=SplitEandP(HadronDefinition, string, newString);

       delete newString; newString=0;
	
       G4KineticTrack * Hadron =0;
       if ( HadronMomentum != 0 ) {

#ifdef debug_VStringDecay                     
  G4cout<<"The attempt was successful"<<G4endl;
#endif
	   G4ThreeVector   Pos;
	   Hadron = new G4KineticTrack(HadronDefinition, 0,Pos, *HadronMomentum);

           if(HadronDefinition->GetPDGIsospin() >= 1.5) {                 // Uzhi June 2016
             G4LorentzVector* Tmp = new G4LorentzVector(*HadronMomentum);
             G4double HadMass2 = Tmp->mag2();

	     Tmp->setPx(0.); Tmp->setPy(0.);
             Tmp->setE(std::sqrt(Tmp->vect().mag2()+HadMass2)); 
	     newString=new G4FragmentingString(*string,newStringEnd,
	   				Tmp);
	     delete Tmp;                                                  // Uzhi June 2016
           } else {
	     newString=new G4FragmentingString(*string,newStringEnd,
	   				HadronMomentum);
           }
	   delete HadronMomentum;
       }
       else
       {

#ifdef debug_VStringDecay
  G4cout<<"The attempt was not successful !!!"<<G4endl;
#endif
       }

#ifdef debug_VStringDecay
  G4cout<<"End SplitUP (G4VLongitudinalStringDecay) ====================="<<G4endl;
#endif

       return Hadron;
}
*/                                                        // Uzhi 28 June 2016, shifted to FTF and QGSM
//--------------------------------------------------------------------------------------

G4ParticleDefinition *
		G4VLongitudinalStringDecay::QuarkSplitup(G4ParticleDefinition*
		decay, G4ParticleDefinition *&created)
{
    G4int IsParticle=(decay->GetPDGEncoding()>0) ? -1 : +1; // if we have a quark, 
                                                            // we need antiquark 
                                                            // (or diquark)
    pDefPair QuarkPair = CreatePartonPair(IsParticle);
    created = QuarkPair.second;
    return hadronizer->Build(QuarkPair.first, decay);

}

//-----------------------------------------------------------------------------

G4int G4VLongitudinalStringDecay::SampleQuarkFlavor(void)
   {
   return (1 + (int)(G4UniformRand()/StrangeSuppress));
   }

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------
G4ThreeVector G4VLongitudinalStringDecay::SampleQuarkPt(G4double ptMax)
   {
   G4double Pt;
   if ( ptMax < 0 ) {
      // sample full gaussian
      Pt = -G4Log(G4UniformRand());
   } else {
      // sample in limited range
      Pt = -G4Log(G4RandFlat::shoot(G4Exp(-sqr(ptMax)/sqr(SigmaQT)), 1.));
   }
   Pt = SigmaQT * std::sqrt(Pt);
   G4double phi = 2.*pi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
   }

//******************************************************************************

void G4VLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector* Hadrons)
   {

//   `yo-yo` formation time
//   const G4double kappa = 1.0 * GeV/fermi/4.;      
   G4double kappa = GetStringTensionParameter();
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
      Hadrons->operator[](c1)->SetFormationTime(
(theInitialStringMass - 2.*SumPz + HadronE - HadronPz)/(2.*kappa)/c_light); 

      G4ThreeVector aPosition(0, 0,     
(theInitialStringMass - 2.*SumE  - HadronE + HadronPz)/(2.*kappa));
      Hadrons->operator[](c1)->SetPosition(aPosition);

      } 
   }

//-----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------

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
		throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetVectorMesonProbability after FragmentString() not allowed");
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
		throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability after FragmentString() not allowed");
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
		throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetScalarMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetScalarMesonMixings( argument Vector too small");
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
		throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetVectorMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "G4VLongitudinalStringDecay::SetVectorMesonMixings( argument Vector too small");
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

//-------------------------------------------------------------------------------------------
void G4VLongitudinalStringDecay::SetStringTensionParameter(G4double aValue)
{
          Kappa = aValue * GeV/fermi;
}	
//**************************************************************************************

