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

#include "G4HadronicException.hh" 

//------------------------debug switches
//#define debug_VStringDecay
//#define debug_heavyHadrons

//******************************************************************************
// Constructors

G4VLongitudinalStringDecay::G4VLongitudinalStringDecay(const G4String& name) 
  : G4HadronicInteraction(name), ProbCCbar(0.0), ProbBBbar(0.0)
{
   MassCut = 210.0*MeV;   // Mpi + Delta

   StringLoopInterrupt  = 1000;
   ClusterLoopInterrupt =  500;

   // Changable Parameters below.
   SigmaQT = 0.5 * GeV;
   
   StrangeSuppress  = 0.44;    // =0.27/2.27 suppression of strange quark pair production, ie. u:d:s=1:1:0.27
   DiquarkSuppress  = 0.07;    // Probability of qq-qqbar pair production
   DiquarkBreakProb = 0.1;     // Probability of (qq)->h+(qq)'
   
   //... pspin_meson is probability to create pseudo-scalar meson 
   pspin_meson.resize(3);
   pspin_meson[0] = 0.5;  // u or d + anti-u or anti-d
   pspin_meson[1] = 0.4;  // one of the quark is strange, or charm, or bottom
   pspin_meson[2] = 0.3;  // both of the quark are strange, or charm, or bottom
   
   //... pspin_barion is probability to create 1/2 barion 
   pspin_barion = 0.5;

   //... vectorMesonMix[] is quark mixing parameters for vector mesons (Variable spin = 3)
   vectorMesonMix.resize(6);
   vectorMesonMix[0] = 0.0;
   vectorMesonMix[1] = 0.5;
   vectorMesonMix[2] = 0.0;
   vectorMesonMix[3] = 0.5;
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

   SetProbCCbar(0.0);  // Probability of CCbar pair creation
   SetProbEta_c(0.1);  // Mixing of Eta_c and J/Psi
   SetProbBBbar(0.0);  // Probability of BBbar pair creation
   SetProbEta_b(0.0);  // Mixing of Eta_b and Upsilon_b

   // Parameters may be changed until the first fragmentation starts
   PastInitPhase=false;
   hadronizer = new G4HadronBuilder( pspin_meson, pspin_barion, scalarMesonMix, vectorMesonMix,
                                     ProbEta_c, ProbEta_b );

   MaxMass=-350.0*GeV;  // If there will be a particle with mass larger than Higgs the value must be changed.

   SetMinMasses();  // Re-calculation of minimal mass of strings and weights of particles in 2-part. decays

   Kappa = 1.0 * GeV/fermi;
   DecayQuark = NewQuark = 0;
}

G4VLongitudinalStringDecay::~G4VLongitudinalStringDecay()
{
   delete hadronizer;
}

G4HadFinalState* 
G4VLongitudinalStringDecay::ApplyYourself(const G4HadProjectile&, G4Nucleus&)
{
  return nullptr;
}

//=============================================================================

// For changing Mass Cut used for selection of very small mass strings
void G4VLongitudinalStringDecay::SetMassCut(G4double aValue){ MassCut=aValue; }
G4double G4VLongitudinalStringDecay::GetMassCut() { return MassCut; }

//-----------------------------------------------------------------------------

// For handling a string with very low mass

G4KineticTrackVector* G4VLongitudinalStringDecay::ProduceOneHadron(const G4ExcitedString * const string)
{
        G4KineticTrackVector* result = nullptr; 
        pDefPair hadrons( nullptr, nullptr );
        G4FragmentingString aString( *string );

        #ifdef debug_VStringDecay
        G4cout<<"G4VLongitudinalStringDecay::ProduceOneHadron: PossibleHmass StrMass "
              <<aString.Mass()<<" MassCut "<<MassCut<<G4endl;
        #endif
        
        SetMinimalStringMass( &aString );
        PossibleHadronMass( &aString, 0, &hadrons );
        result = new G4KineticTrackVector;
        if ( hadrons.first != nullptr ) {       
           if ( hadrons.second == nullptr ) {
               // Substitute string by light hadron, Note that Energy is not conserved here!

               #ifdef debug_VStringDecay
               G4cout << "VlongSD Warning replacing string by single hadron (G4VLongitudinalStringDecay)" <<G4endl;
               G4cout << hadrons.first->GetParticleName()<<G4endl
                      << "string .. " << string->Get4Momentum() << " " 
                      << string->Get4Momentum().m() << G4endl;
               #endif           

               G4ThreeVector   Mom3 = string->Get4Momentum().vect();
               G4LorentzVector Mom( Mom3, std::sqrt( Mom3.mag2() + sqr( hadrons.first->GetPDGMass() ) ) );
               result->push_back( new G4KineticTrack( hadrons.first, 0, string->GetPosition(), Mom ) );
           } else {
               //... string was qq--qqbar type: Build two stable hadrons,

               #ifdef debug_VStringDecay
               G4cout << "VlongSD Warning replacing qq-qqbar string by TWO hadrons (G4VLongitudinalStringDecay)" 
                      << hadrons.first->GetParticleName() << " / " 
                      << hadrons.second->GetParticleName()
                      << "string .. " << string->Get4Momentum() << " " 
                      << string->Get4Momentum().m() << G4endl;
               #endif

               G4LorentzVector  Mom1, Mom2;
               Sample4Momentum( &Mom1, hadrons.first->GetPDGMass(), 
                                &Mom2, hadrons.second->GetPDGMass(),
                                string->Get4Momentum().mag() );

               result->push_back( new G4KineticTrack( hadrons.first,  0, string->GetPosition(), Mom1 ) );
               result->push_back( new G4KineticTrack( hadrons.second, 0, string->GetPosition(), Mom2 ) );

               G4ThreeVector Velocity = string->Get4Momentum().boostVector();
               result->Boost(Velocity);          
           }
        }
        return result;
}

//----------------------------------------------------------------------------------------

G4double G4VLongitudinalStringDecay::PossibleHadronMass( const G4FragmentingString * const string,
                                                         Pcreate build, pDefPair * pdefs )
{
        G4double mass = 0.0;

	if ( build==0 ) build=&G4HadronBuilder::BuildLowSpin;

        G4ParticleDefinition* Hadron1 = nullptr;
	G4ParticleDefinition* Hadron2 = nullptr;

        if (!string->IsAFourQuarkString() )
        {
           // spin 0 meson or spin 1/2 barion will be built

           Hadron1 = (hadronizer->*build)(string->GetLeftParton(), string->GetRightParton());
           #ifdef debug_VStringDecay
	   G4cout<<"VlongSD PossibleHadronMass"<<G4endl;
           G4cout<<"VlongSD Quarks at the string ends "<<string->GetLeftParton()->GetParticleName()
                 <<" "<<string->GetRightParton()->GetParticleName()<<G4endl;
           if ( Hadron1 != nullptr) {
             G4cout<<"(G4VLongitudinalStringDecay) Hadron "<<Hadron1->GetParticleName()
                   <<" "<<Hadron1->GetPDGMass()<<G4endl;
           }
           #endif
           if ( Hadron1 != nullptr) { mass = (Hadron1)->GetPDGMass();}
           else                  { mass = MaxMass;}
        } else
        {
           //... string is qq--qqbar: Build two stable hadrons,

           #ifdef debug_VStringDecay
           G4cout<<"VlongSD PossibleHadronMass"<<G4endl;
           G4cout<<"VlongSD string is qq--qqbar: Build two stable hadrons"<<G4endl; 
           #endif

           G4double StringMass   = string->Mass();
           G4int cClusterInterrupt = 0;
           do
           {
             if (cClusterInterrupt++ >= ClusterLoopInterrupt) return false;

             G4int LeftQuark1= string->GetLeftParton()->GetPDGEncoding()/1000;
             G4int LeftQuark2=(string->GetLeftParton()->GetPDGEncoding()/100)%10;

             G4int RightQuark1= string->GetRightParton()->GetPDGEncoding()/1000;
             G4int RightQuark2=(string->GetRightParton()->GetPDGEncoding()/100)%10;

             if (G4UniformRand()<0.5) {
               Hadron1 =hadronizer->Build(FindParticle(LeftQuark1), FindParticle(RightQuark1));
               Hadron2 =hadronizer->Build(FindParticle(LeftQuark2), FindParticle(RightQuark2));
             } else {
               Hadron1 =hadronizer->Build(FindParticle(LeftQuark1), FindParticle(RightQuark2));
               Hadron2 =hadronizer->Build(FindParticle(LeftQuark2), FindParticle(RightQuark1));
             }
             //... repeat procedure, if mass of cluster is too low to produce hadrons
             //... ClusterMassCut = 0.15*GeV model parameter
           }
           while ( Hadron1 == nullptr || Hadron2 == nullptr ||
                   ( StringMass <= Hadron1->GetPDGMass() + Hadron2->GetPDGMass() ) );

	   mass = (Hadron1)->GetPDGMass() + (Hadron2)->GetPDGMass();
        }

        #ifdef debug_VStringDecay
        G4cout<<"VlongSD *Hadrons 1 and 2, proposed mass "<<Hadron1<<" "<<Hadron2<<" "<<mass<<G4endl; 
        #endif
	
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
  /*
  G4cout<<Encoding<<" G4VLongitudinalStringDecay::FindParticle Check di-quarks *******************"<<G4endl;
  for (G4int i=4; i<6;i++){
    for (G4int j=1;j<6;j++){
      G4cout<<i<<" "<<j<<" ";
      G4int Code = 1000 * i + 100 * j +1;
      G4ParticleDefinition* ptr1 = G4ParticleTable::GetParticleTable()->FindParticle(Code);
      Code +=2;
      G4ParticleDefinition* ptr2 = G4ParticleTable::GetParticleTable()->FindParticle(Code);
      G4cout<<"Code "<<Code - 2<<" ptr "<<ptr1<<" :: Code "<<Code<<" ptr "<<ptr2<<G4endl;
    }
    G4cout<<G4endl;
  }
  */

  G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);

  if (ptr == nullptr)
  {
     for (size_t i=0; i < NewParticles.size(); i++)
     {
       if ( Encoding == NewParticles[i]->GetPDGEncoding() ) { ptr = NewParticles[i]; return ptr;}
     }
  }

  return ptr;    
}

//*********************************************************************************
//   For decision on continue or stop string fragmentation
//   virtual G4bool StopFragmenting(const G4FragmentingString  * const string)=0;
//   virtual G4bool IsItFragmentable(const G4FragmentingString * const string)=0;
//
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

G4ParticleDefinition * G4VLongitudinalStringDecay::QuarkSplitup( G4ParticleDefinition* decay,
                                                                 G4ParticleDefinition *&created )
{
   #ifdef debug_VStringDecay
   G4cout<<"VlongSD QuarkSplitup: quark ID "<<decay->GetPDGEncoding()<<G4endl; 
   #endif
   
   G4int IsParticle=(decay->GetPDGEncoding()>0) ? -1 : +1;  // if we have a quark, we need antiquark (or diquark)

   pDefPair QuarkPair = CreatePartonPair(IsParticle);
   created = QuarkPair.second;

   DecayQuark = decay->GetPDGEncoding();
   NewQuark   = created->GetPDGEncoding();

   #ifdef debug_VStringDecay
   G4cout<<"VlongSD QuarkSplitup: "<<decay->GetPDGEncoding()<<" -> "<<QuarkPair.second->GetPDGEncoding()<<G4endl;
   G4cout<<"hadronizer->Build(QuarkPair.first, decay)"<<G4endl;
   #endif
   
   return hadronizer->Build(QuarkPair.first, decay);
}

//-----------------------------------------------------------------------------

G4VLongitudinalStringDecay::pDefPair G4VLongitudinalStringDecay::
CreatePartonPair(G4int NeedParticle,G4bool AllowDiquarks)
{
    //  NeedParticle = +1 for Particle, -1 for Antiparticle
    if ( AllowDiquarks && G4UniformRand() < DiquarkSuppress )
    {
      // Create a Diquark - AntiDiquark pair , first in pair is anti to IsParticle
      #ifdef debug_VStringDecay
      G4cout<<"VlongSD Create a Diquark - AntiDiquark pair"<<G4endl;
      #endif
      G4int q1(0), q2(0), spin(0), PDGcode(0);

      q1  = SampleQuarkFlavor();
      q2  = SampleQuarkFlavor();

      spin = (q1 != q2 && G4UniformRand() <= 0.5)? 1 : 3;
                                     //   convention: quark with higher PDG number is first
      PDGcode = (std::max(q1,q2) * 1000 + std::min(q1,q2) * 100 + spin) * NeedParticle;

      return pDefPair (FindParticle(-PDGcode),FindParticle(PDGcode));

    } else {
      // Create a Quark - AntiQuark pair, first in pair  IsParticle
      #ifdef debug_VStringDecay
      G4cout<<"VlongSD Create a Quark - AntiQuark pair"<<G4endl; 
      #endif
      G4int PDGcode=SampleQuarkFlavor()*NeedParticle;
      return pDefPair (FindParticle(PDGcode),FindParticle(-PDGcode));
    }
}

//-----------------------------------------------------------------------------

G4int G4VLongitudinalStringDecay::SampleQuarkFlavor(void)
{
   G4int  quark(1);
   G4double ksi = G4UniformRand();
   if ( ksi < ProbCB ) {
      if ( ksi < ProbCCbar ) {quark = 4;}   // c quark
      else                   {quark = 5;}   // b quark
      #ifdef debug_heavyHadrons
      G4cout << "G4VLongitudinalStringDecay::SampleQuarkFlavor : sampled from the vacuum HEAVY quark = "
	     << quark << G4endl;
      #endif
   } else {
     quark = 1 + (int)(G4UniformRand()/StrangeSuppress);
   }
   #ifdef debug_VStringDecay
   G4cout<<"VlongSD SampleQuarkFlavor "<<quark<<" (ProbCB ProbCCbar ProbBBbar "<<ProbCB
         <<" "<<ProbCCbar<<" "<<ProbBBbar<<" )"<<G4endl; 
   #endif
   return quark;
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
      G4double q = ptMax/SigmaQT;
      G4double ymin = (q > 20.) ? 0.0 : G4Exp(-q*q); 
      Pt = -G4Log(G4RandFlat::shoot(ymin, 1.));
   }
   Pt = SigmaQT * std::sqrt(Pt);
   G4double phi = 2.*pi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
}

//******************************************************************************

void G4VLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, 
                                                             G4KineticTrackVector* Hadrons)
{
   //   `yo-yo` formation time
   //   const G4double kappa = 1.0 * GeV/fermi/4.;      
   G4double kappa = GetStringTensionParameter();
   for (size_t c1 = 0; c1 < Hadrons->size(); c1++)
   {
      G4double SumPz = 0; 
      G4double SumE  = 0;
      for (size_t c2 = 0; c2 < c1; c2++)
      {
         SumPz += Hadrons->operator[](c2)->Get4Momentum().pz();
         SumE  += Hadrons->operator[](c2)->Get4Momentum().e();   
      } 
      G4double HadronE  = Hadrons->operator[](c1)->Get4Momentum().e();
      G4double HadronPz = Hadrons->operator[](c1)->Get4Momentum().pz();
      Hadrons->operator[](c1)->SetFormationTime(
        (theInitialStringMass - 2.*SumPz + HadronE - HadronPz ) / (2.*kappa) / c_light ); 
      G4ThreeVector aPosition( 0, 0,
        (theInitialStringMass - 2.*SumE  - HadronE + HadronPz) / (2.*kappa) );
      Hadrons->operator[](c1)->SetPosition(aPosition);
   }
}

//-----------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetSigmaTransverseMomentum(G4double aValue)
{
   if ( PastInitPhase ) {
     throw G4HadronicException(__FILE__, __LINE__, 
       "G4VLongitudinalStringDecay::SetSigmaTransverseMomentum after FragmentString() not allowed");
   } else {
     SigmaQT = aValue;
   }
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetStrangenessSuppression(G4double aValue)
{
   StrangeSuppress = aValue;
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetDiquarkSuppression(G4double aValue)
{
   DiquarkSuppress = aValue;
}

//----------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetDiquarkBreakProbability(G4double aValue)
{
  if ( PastInitPhase ) {
    throw G4HadronicException(__FILE__, __LINE__, 
      "G4VLongitudinalStringDecay::SetDiquarkBreakProbability after FragmentString() not allowed");
  } else {
    DiquarkBreakProb = aValue;
  }
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability(G4double aValue)
{
  if ( PastInitPhase ) {
    throw G4HadronicException(__FILE__, __LINE__, 
      "G4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability after FragmentString() not allowed");
  } else {
    pspin_barion = aValue;
    delete hadronizer;
    hadronizer = new G4HadronBuilder( pspin_meson, pspin_barion, scalarMesonMix, vectorMesonMix, 
                                      ProbEta_c, ProbEta_b );
  }
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetScalarMesonMixings(std::vector<G4double> aVector)
{
  if ( PastInitPhase ) {
    throw G4HadronicException(__FILE__, __LINE__, 
      "G4VLongitudinalStringDecay::SetScalarMesonMixings after FragmentString() not allowed");
  } else {
    if ( aVector.size() < 6 ) 
      throw G4HadronicException(__FILE__, __LINE__, 
        "G4VLongitudinalStringDecay::SetScalarMesonMixings( argument Vector too small");
    scalarMesonMix[0] = aVector[0];
    scalarMesonMix[1] = aVector[1];
    scalarMesonMix[2] = aVector[2];
    scalarMesonMix[3] = aVector[3];
    scalarMesonMix[4] = aVector[4];
    scalarMesonMix[5] = aVector[5];
    delete hadronizer;
    hadronizer = new G4HadronBuilder( pspin_meson, pspin_barion, scalarMesonMix, vectorMesonMix, 
                                      ProbEta_c, ProbEta_b );
  }
}

//----------------------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetVectorMesonMixings(std::vector<G4double> aVector)
{
  if ( PastInitPhase ) {
    throw G4HadronicException(__FILE__, __LINE__, 
      "G4VLongitudinalStringDecay::SetVectorMesonMixings after FragmentString() not allowed");
  } else {
    if ( aVector.size() < 6 ) 
      throw G4HadronicException(__FILE__, __LINE__, 
        "G4VLongitudinalStringDecay::SetVectorMesonMixings( argument Vector too small");
    vectorMesonMix[0] = aVector[0];
    vectorMesonMix[1] = aVector[1];
    vectorMesonMix[2] = aVector[2];
    vectorMesonMix[3] = aVector[3];
    vectorMesonMix[4] = aVector[4];
    vectorMesonMix[5] = aVector[5];
    delete hadronizer;
    hadronizer = new G4HadronBuilder( pspin_meson, pspin_barion, scalarMesonMix, vectorMesonMix, 
                                      ProbEta_c, ProbEta_b );
  }
}

//-------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetProbCCbar(G4double aValue)
{
   ProbCCbar = aValue;
   ProbCB = ProbCCbar + ProbBBbar;
}

//-------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetProbEta_c(G4double aValue)
{
   ProbEta_c = aValue;
}

//-------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetProbBBbar(G4double aValue)
{
   ProbBBbar = aValue;
   ProbCB = ProbCCbar + ProbBBbar;
}

//-------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetProbEta_b(G4double aValue)
{
   ProbEta_b = aValue;
}

//-------------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetStringTensionParameter(G4double aValue)
{
   Kappa = aValue * GeV/fermi;
}	

//-----------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetMinMasses()
{
    // ------ For estimation of a minimal string mass ---------------
    Mass_of_light_quark =140.*MeV;
    Mass_of_s_quark     =500.*MeV;
    Mass_of_c_quark     =1600.*MeV;
    Mass_of_b_quark     =4500.*MeV;
    Mass_of_string_junction=720.*MeV;

    // ---------------- Determination of minimal mass of q-qbar strings -------------------
    G4ParticleDefinition * hadron1;    G4int Code1;
    G4ParticleDefinition * hadron2;    G4int Code2;
    for (G4int i=1; i < 6; i++) {
        Code1 = 100*i + 10*1 + 1;
        hadron1 = FindParticle(Code1);

        if (hadron1 != nullptr) {
           for (G4int j=1; j < 6; j++) {
               Code2 = 100*j + 10*1 + 1;
               hadron2 = FindParticle(Code2);
               if (hadron2 != nullptr) {
                 minMassQQbarStr[i-1][j-1] = hadron1->GetPDGMass() + hadron2->GetPDGMass() + 70.0 * MeV;
               }
           } 
        }
    }

    minMassQQbarStr[1][1] = minMassQQbarStr[0][0];   // u-ubar = 0.5 Pi0 + 0.24 Eta + 0.25 Eta'

    // ---------------- Determination of minimal mass of qq-q strings -------------------
    G4ParticleDefinition * hadron3;
    G4int kfla, kflb;
    //  MaxMass = -350.0*GeV;   // If there will be a particle with mass larger than Higgs the value must be changed.

    for (G4int i=1; i < 6; i++) {   //i=1
        Code1 = 100*i + 10*1 + 1;
        hadron1 = FindParticle(Code1);
        for (G4int j=1; j < 6; j++) {
            for (G4int k=1; k < 6; k++) {
                kfla = std::max(j,k);
                kflb = std::min(j,k);

		// Add d-quark
                Code2 = 1000*kfla + 100*kflb + 10*1 + 2;
		if ( (j == 1) && (k==1)) Code2 = 1000*2 + 100*1 + 10*1 + 2; // In the case - add u-quark.

                hadron2 = G4ParticleTable::GetParticleTable()->FindParticle(Code2);
                hadron3 = G4ParticleTable::GetParticleTable()->FindParticle(Code2 + 2);

                if ((hadron2 == nullptr) && (hadron3 == nullptr)) {minMassQDiQStr[i-1][j-1][k-1] = MaxMass; continue;};

                if ((hadron2 != nullptr) && (hadron3 != nullptr)) {
                   if (hadron2->GetPDGMass() > hadron3->GetPDGMass() ) { hadron2 = hadron3; }
                };

                if ((hadron2 != nullptr) && (hadron3 == nullptr)) {};

                if ((hadron2 == nullptr) && (hadron3 != nullptr)) {hadron2 = hadron3;};

                minMassQDiQStr[i-1][j-1][k-1] = hadron1->GetPDGMass() + hadron2->GetPDGMass() + 70.0 * MeV;
            }
        }
    }

    // ------ An estimated minimal string mass ----------------------
    MinimalStringMass  = 0.;
    MinimalStringMass2 = 0.;
    // q charges  d               u                s               c                b
    Qcharge[0] = -1; Qcharge[1] = 2; Qcharge[2] = -1; Qcharge[3] = 2; Qcharge[4] = -1;

    // For treating of small string decays
    for (G4int i=0; i<5; i++)
    {  for (G4int j=0; j<5; j++)
       {  for (G4int k=0; k<7; k++)
          {
            Meson[i][j][k]=0; MesonWeight[i][j][k]=0.;
          }
       }
    }
    //--------------------------
    G4int StrangeQ = 0;
    G4int StrangeAQ = 0;
    for (G4int i=0; i<5; i++)
    {
       if( i >= 2 ) StrangeQ=1;
       for (G4int j=0; j<5; j++)
       { 
         StrangeAQ = 0;
         if( j >= 2 ) StrangeAQ=1;
         Meson[i][j][0]       = 100 * (std::max(i,j)+1) + 10 * (std::min(i,j)+1) + 1; // Scalar meson
         MesonWeight[i][j][0] = (   pspin_meson[StrangeQ + StrangeAQ]);
         Meson[i][j][1]       = 100 * (std::max(i,j)+1) + 10 * (std::min(i,j)+1) + 3; // Vector meson
         MesonWeight[i][j][1] = (1.-pspin_meson[StrangeQ + StrangeAQ]);
       }
    }

    //qqs                                                                                                         indexes
    //dd1 -> scalarMesonMix[0] * 111 + (1-scalarMesonMix[0]-scalarMesonMix[1]) * 221 + scalarMesonMix[1] * 331     (000)
    //dd1 ->                     Pi0                                             Eta                       Eta'

    Meson[0][0][0] = 111; MesonWeight[0][0][0] = (   pspin_meson[0]) * (  scalarMesonMix[0]                  );  // Pi0
    Meson[0][0][2] = 221; MesonWeight[0][0][3] = (   pspin_meson[0]) * (1-scalarMesonMix[0]-scalarMesonMix[1]);  // Eta
    Meson[0][0][3] = 331; MesonWeight[0][0][4] = (   pspin_meson[0]) * (                    scalarMesonMix[1]);  // Eta'

    //dd3 -> (1-vectorMesonMix[1] * 113 + vectorMesonMix[1] * 223                                                  (001)
    //dd3 ->                       rho_0                     omega

    Meson[0][0][1] = 113; MesonWeight[0][0][1] = (1.-pspin_meson[0]) * (1-vectorMesonMix[1]);                    // Rho
    Meson[0][0][4] = 223; MesonWeight[0][0][4] = (1.-pspin_meson[0]) * (  vectorMesonMix[1]);                    // omega

    //uu1 -> scalarMesonMix[0] * 111 + (1-scalarMesonMix[0]-scalarMesonMix[1]) * 221 + scalarMesonMix[1] * 331     (110)
    //uu1 ->                     Pi0                                             Eta                       Eta'

    Meson[1][1][0] = 111; MesonWeight[1][1][0] = (   pspin_meson[0]) * (  scalarMesonMix[0]                  );  // Pi0
    Meson[1][1][2] = 221; MesonWeight[1][1][2] = (   pspin_meson[0]) * (1-scalarMesonMix[0]-scalarMesonMix[1]);  // Eta
    Meson[1][1][3] = 331; MesonWeight[1][1][3] = (   pspin_meson[0]) * (                    scalarMesonMix[1]);  // Eta'

    //uu3 -> (1-vectorMesonMix[1]) * 113 + vectorMesonMix[1] * 223                                                 (111)
    //uu3 ->                        rho_0                     omega

    Meson[1][1][1] = 113; MesonWeight[1][1][1] = (1.-pspin_meson[0]) * (1-vectorMesonMix[1]);                    // Rho
    Meson[1][1][4] = 223; MesonWeight[1][1][4] = (1.-pspin_meson[0]) * (  vectorMesonMix[1]);                    // omega

    //ss1     ->                                             (1-scalarMesonMix[5]) * 221 + scalarMesonMix[5] * 331   (220)
    //ss1     ->                                                                     Eta                       Eta'

    Meson[2][2][0] = 221; MesonWeight[2][2][0] = (   pspin_meson[2]) * (1-scalarMesonMix[5]                  );  // Eta
    Meson[2][2][2] = 331; MesonWeight[2][2][2] = (   pspin_meson[2]) * (                    scalarMesonMix[5]);  // Eta'

    //ss3     ->                                                                           vectorMesonMix[5] * 333   (221)
    //ss3     ->                                                                                               phi

    Meson[2][2][1] = 333; MesonWeight[2][2][1] = (1.-pspin_meson[2]) * (                 vectorMesonMix[5]);  // phi

    //cc1     ->    ProbEta_c /(1-pspin_meson) 441  (330) Probability of Eta_c
    //cc3     -> (1-ProbEta_c)/(  pspin_meson) 443  (331) Probability of J/Psi

    //bb1     ->    ProbEta_b /pspin_meson 551  (440) Probability of Eta_b
    //bb3     -> (1-ProbEta_b)/pspin_meson 553  (441) Probability of Upsilon

    if ( pspin_meson[2] != 0. ) {
       Meson[3][3][0] *= (    ProbEta_c)/(   pspin_meson[2]);   // Eta_c
       Meson[3][3][1] *= (1.0-ProbEta_c)/(1.-pspin_meson[2]);   // J/Psi

       Meson[4][4][0] *= (    ProbEta_b)/(   pspin_meson[2]);   // Eta_b
       Meson[4][4][1] *= (1.0-ProbEta_b)/(1.-pspin_meson[2]);   // Upsilon
    }

    //--------------------------

    for (G4int i=0; i<5; i++)
    {  for (G4int j=0; j<5; j++)
       {  for (G4int k=0; k<5; k++)
          {  for (G4int l=0; l<4; l++)
             { Baryon[i][j][k][l]=0; BaryonWeight[i][j][k][l]=0.;}
          }
       }
    }

          kfla =0;  kflb =0;
    G4int                   kflc(0), kfld(0), kfle(0), kflf(0);
    for (G4int i=0; i<5; i++)
    {  for (G4int j=0; j<5; j++)
       {  for (G4int k=0; k<5; k++)
          {  
           kfla = i+1; kflb = j+1; kflc = k+1;
	   kfld = std::max(kfla,kflb);
	   kfld = std::max(kfld,kflc);

	   kflf = std::min(kfla,kflb);
	   kflf = std::min(kflf,kflc);

           kfle = kfla + kflb + kflc - kfld - kflf;

           Baryon[i][j][k][0]       = 1000 * kfld + 100 * kfle + 10 * kflf + 2; // spin=1/2
           BaryonWeight[i][j][k][0] = (   pspin_barion);
           Baryon[i][j][k][1]       = 1000 * kfld + 100 * kfle + 10 * kflf + 4; // spin=3/2
           BaryonWeight[i][j][k][1] = (1.-pspin_barion);
          }
       }
    }

    // Delta-  ddd - only 1114
    Baryon[0][0][0][0] = 1114;    BaryonWeight[0][0][0][0] = 1.0; 
    Baryon[0][0][0][1] =    0;    BaryonWeight[0][0][0][1] = 0.0; 

    // Delta++ uuu - only 2224
    Baryon[1][1][1][0] = 2224;    BaryonWeight[1][1][1][0] = 1.0; 
    Baryon[1][1][1][1] =    0;    BaryonWeight[1][1][1][1] = 0.0; 

    // Omega- sss - only 3334
    Baryon[2][2][2][0] = 3334;    BaryonWeight[2][2][2][0] = 1.0; 
    Baryon[2][2][2][1] =    0;    BaryonWeight[2][2][2][1] = 0.0; 

    // Omega_cc++ ccc - only 4444
    Baryon[3][3][3][0] = 4444;    BaryonWeight[3][3][3][0] = 1.0; 
    Baryon[3][3][3][1] =    0;    BaryonWeight[3][3][3][1] = 0.0; 

    // Omega_bb-  bbb - only 5554
    Baryon[4][4][4][0] = 5554;    BaryonWeight[4][4][4][0] = 1.0; 
    Baryon[4][4][4][1] =    0;    BaryonWeight[4][4][4][1] = 0.0; 

    // Lambda/Sigma0 sud - 3122/3212
    Baryon[0][1][2][0] = 3122;    BaryonWeight[0][1][2][0] *= 0.5;                  // Lambda
    Baryon[0][2][1][0] = 3122;    BaryonWeight[0][2][1][0] *= 0.5;
    Baryon[1][0][2][0] = 3122;    BaryonWeight[1][0][2][0] *= 0.5;
    Baryon[1][2][0][0] = 3122;    BaryonWeight[1][2][0][0] *= 0.5;
    Baryon[2][0][1][0] = 3122;    BaryonWeight[2][0][1][0] *= 0.5;
    Baryon[2][1][0][0] = 3122;    BaryonWeight[2][1][0][0] *= 0.5;

    Baryon[0][1][2][2] = 3212;    BaryonWeight[0][1][2][2]  = 0.5 * pspin_barion;   // Sigma0
    Baryon[0][2][1][2] = 3212;    BaryonWeight[0][2][1][2]  = 0.5 * pspin_barion;
    Baryon[1][0][2][2] = 3212;    BaryonWeight[1][0][2][2]  = 0.5 * pspin_barion;
    Baryon[1][2][0][2] = 3212;    BaryonWeight[1][2][0][2]  = 0.5 * pspin_barion;
    Baryon[2][0][1][2] = 3212;    BaryonWeight[2][0][1][2]  = 0.5 * pspin_barion;
    Baryon[2][1][0][2] = 3212;    BaryonWeight[2][1][0][2]  = 0.5 * pspin_barion;

    // Lambda_c+/Sigma_c+ cud - 4122/4212
    Baryon[0][1][3][0] = 4122;    BaryonWeight[0][1][3][0] *= 0.5;                  // Lambda_c+
    Baryon[0][3][1][0] = 4122;    BaryonWeight[0][3][1][0] *= 0.5;
    Baryon[1][0][3][0] = 4122;    BaryonWeight[1][0][3][0] *= 0.5;
    Baryon[1][3][0][0] = 4122;    BaryonWeight[1][3][0][0] *= 0.5;
    Baryon[3][0][1][0] = 4122;    BaryonWeight[3][0][1][0] *= 0.5;
    Baryon[3][1][0][0] = 4122;    BaryonWeight[3][1][0][0] *= 0.5;

    Baryon[0][1][3][2] = 4212;    BaryonWeight[0][1][3][2]  = 0.5 * pspin_barion;   // SigmaC+
    Baryon[0][3][1][2] = 4212;    BaryonWeight[0][3][1][2]  = 0.5 * pspin_barion;
    Baryon[1][0][3][2] = 4212;    BaryonWeight[1][0][3][2]  = 0.5 * pspin_barion;
    Baryon[1][3][0][2] = 4212;    BaryonWeight[1][3][0][2]  = 0.5 * pspin_barion;
    Baryon[3][0][1][2] = 4212;    BaryonWeight[3][0][1][2]  = 0.5 * pspin_barion;
    Baryon[3][1][0][2] = 4212;    BaryonWeight[3][1][0][2]  = 0.5 * pspin_barion;

    // Xi_c+/Xi_c+' cus - 4232/4322
    Baryon[1][2][3][0] = 4232;    BaryonWeight[1][2][3][0] *= 0.5;                  // Xi_c+
    Baryon[1][3][2][0] = 4232;    BaryonWeight[1][3][2][0] *= 0.5;
    Baryon[2][1][3][0] = 4232;    BaryonWeight[2][1][3][0] *= 0.5;
    Baryon[2][3][1][0] = 4232;    BaryonWeight[2][3][1][0] *= 0.5;
    Baryon[3][1][2][0] = 4232;    BaryonWeight[3][1][2][0] *= 0.5;
    Baryon[3][2][1][0] = 4232;    BaryonWeight[3][2][1][0] *= 0.5;

    Baryon[1][2][3][2] = 4322;    BaryonWeight[1][2][3][2]  = 0.5 * pspin_barion;   // Xi_c+'
    Baryon[1][3][2][2] = 4322;    BaryonWeight[1][3][2][2]  = 0.5 * pspin_barion;
    Baryon[2][1][3][2] = 4322;    BaryonWeight[2][1][3][2]  = 0.5 * pspin_barion;
    Baryon[2][3][1][2] = 4322;    BaryonWeight[2][3][1][2]  = 0.5 * pspin_barion;
    Baryon[3][1][2][2] = 4322;    BaryonWeight[3][1][2][2]  = 0.5 * pspin_barion;
    Baryon[3][2][1][2] = 4322;    BaryonWeight[3][2][1][2]  = 0.5 * pspin_barion;

    // Xi_c0/Xi_c0' cus - 4132/4312
    Baryon[0][2][3][0] = 4132;    BaryonWeight[0][2][3][0] *= 0.5;                  // Xi_c0
    Baryon[0][3][2][0] = 4132;    BaryonWeight[0][3][2][0] *= 0.5;
    Baryon[2][0][3][0] = 4132;    BaryonWeight[2][0][3][0] *= 0.5;
    Baryon[2][3][0][0] = 4132;    BaryonWeight[2][3][0][0] *= 0.5;
    Baryon[3][0][2][0] = 4132;    BaryonWeight[3][0][2][0] *= 0.5;
    Baryon[3][2][0][0] = 4132;    BaryonWeight[3][2][0][0] *= 0.5;

    Baryon[0][2][3][2] = 4312;    BaryonWeight[0][2][3][2]  = 0.5 * pspin_barion;   // Xi_c0'
    Baryon[0][3][2][2] = 4312;    BaryonWeight[0][3][2][2]  = 0.5 * pspin_barion;
    Baryon[2][0][3][2] = 4312;    BaryonWeight[2][0][3][2]  = 0.5 * pspin_barion;
    Baryon[2][3][0][2] = 4312;    BaryonWeight[2][3][0][2]  = 0.5 * pspin_barion;
    Baryon[3][0][2][2] = 4312;    BaryonWeight[3][0][2][2]  = 0.5 * pspin_barion;
    Baryon[3][2][0][2] = 4312;    BaryonWeight[3][2][0][2]  = 0.5 * pspin_barion;

    // Lambda_b0/Sigma_b0 bud - 5122/5212
    Baryon[0][1][4][0] = 5122;    BaryonWeight[0][1][4][0] *= 0.5;                  // Lambda_b0
    Baryon[0][4][1][0] = 5122;    BaryonWeight[0][4][1][0] *= 0.5;
    Baryon[1][0][4][0] = 5122;    BaryonWeight[1][0][4][0] *= 0.5;
    Baryon[1][4][0][0] = 5122;    BaryonWeight[1][4][0][0] *= 0.5;
    Baryon[4][0][1][0] = 5122;    BaryonWeight[4][0][1][0] *= 0.5;
    Baryon[4][1][0][0] = 5122;    BaryonWeight[4][1][0][0] *= 0.5;

    Baryon[0][1][4][2] = 5212;    BaryonWeight[0][1][4][2]  = 0.5 * pspin_barion;   // Sigma_b0
    Baryon[0][4][1][2] = 5212;    BaryonWeight[0][4][1][2]  = 0.5 * pspin_barion;
    Baryon[1][0][4][2] = 5212;    BaryonWeight[1][0][4][2]  = 0.5 * pspin_barion;
    Baryon[1][4][0][2] = 5212;    BaryonWeight[1][4][0][2]  = 0.5 * pspin_barion;
    Baryon[4][0][1][2] = 5212;    BaryonWeight[4][0][1][2]  = 0.5 * pspin_barion;
    Baryon[4][1][0][2] = 5212;    BaryonWeight[4][1][0][2]  = 0.5 * pspin_barion;

    // Xi_b0/Xi_b0' bus - 5232/5322
    Baryon[1][2][4][0] = 5232;    BaryonWeight[1][2][4][0] *= 0.5;                  // Xi_b0
    Baryon[1][4][2][0] = 5232;    BaryonWeight[1][4][2][0] *= 0.5;
    Baryon[2][1][4][0] = 5232;    BaryonWeight[2][1][4][0] *= 0.5;
    Baryon[2][4][1][0] = 5232;    BaryonWeight[2][4][1][0] *= 0.5;
    Baryon[4][1][2][0] = 5232;    BaryonWeight[4][1][2][0] *= 0.5;
    Baryon[4][2][1][0] = 5232;    BaryonWeight[4][2][1][0] *= 0.5;

    Baryon[1][2][4][2] = 5322;    BaryonWeight[1][2][4][2]  = 0.5 * pspin_barion;   // Xi_b0'
    Baryon[1][4][2][2] = 5322;    BaryonWeight[1][4][2][2]  = 0.5 * pspin_barion;
    Baryon[2][1][4][2] = 5322;    BaryonWeight[2][1][4][2]  = 0.5 * pspin_barion;
    Baryon[2][4][1][2] = 5322;    BaryonWeight[2][4][1][2]  = 0.5 * pspin_barion;
    Baryon[4][1][2][2] = 5322;    BaryonWeight[4][1][2][2]  = 0.5 * pspin_barion;
    Baryon[4][2][1][2] = 5322;    BaryonWeight[4][2][1][2]  = 0.5 * pspin_barion;

    // Xi_b-/Xi_b-' bus - 5132/5312
    Baryon[0][2][4][0] = 5132;    BaryonWeight[0][2][4][0] *= 0.5;                  // Xi_b-
    Baryon[0][4][2][0] = 5132;    BaryonWeight[0][4][2][0] *= 0.5;
    Baryon[2][0][4][0] = 5132;    BaryonWeight[2][0][4][0] *= 0.5;
    Baryon[2][4][0][0] = 5132;    BaryonWeight[2][4][0][0] *= 0.5;
    Baryon[4][0][2][0] = 5132;    BaryonWeight[4][0][2][0] *= 0.5;
    Baryon[4][2][0][0] = 5132;    BaryonWeight[4][2][0][0] *= 0.5;

    Baryon[0][2][4][2] = 5312;    BaryonWeight[0][2][4][2]  = 0.5 * pspin_barion;   // Xi_b-'
    Baryon[0][4][2][2] = 5312;    BaryonWeight[0][4][2][2]  = 0.5 * pspin_barion;
    Baryon[2][0][4][2] = 5312;    BaryonWeight[2][0][4][2]  = 0.5 * pspin_barion;
    Baryon[2][4][0][2] = 5312;    BaryonWeight[2][4][0][2]  = 0.5 * pspin_barion;
    Baryon[4][0][2][2] = 5312;    BaryonWeight[4][0][2][2]  = 0.5 * pspin_barion;
    Baryon[4][2][0][2] = 5312;    BaryonWeight[4][2][0][2]  = 0.5 * pspin_barion;

    for (G4int i=0; i<5; i++)
    {  for (G4int j=0; j<5; j++)
	  {  for (G4int k=0; k<5; k++)
		 {  for (G4int l=0; l<4; l++)
		    { 
                     G4ParticleDefinition * TestHadron=
                       G4ParticleTable::GetParticleTable()->FindParticle(Baryon[i][j][k][l]);
                     /*
                     G4cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Baryon[i][j][k][l]<<" "<<TestHadron<<" "<<BaryonWeight[i][j][k][l];
                     if (TestHadron != nullptr) G4cout<<" "<<TestHadron->GetParticleName();
                     if ((TestHadron == nullptr)&&(Baryon[i][j][k][l] != 0)) G4cout<<" *****";
                     if ((TestHadron == nullptr)&&(Baryon[i][j][k][l] == 0)) G4cout<<" ---------------";
                     G4cout<<G4endl;
                     */
                     if ((TestHadron == nullptr)&&(Baryon[i][j][k][l] != 0)) Baryon[i][j][k][l] = 0;
                    }
		 }
	  }
    }

    // --------- Probabilities of q-qbar pair productions for kink or gluons.
    G4double ProbUUbar = 0.33;
    Prob_QQbar[0]=ProbUUbar;         // Probability of ddbar production
    Prob_QQbar[1]=ProbUUbar;         // Probability of uubar production
    Prob_QQbar[2]=1.0-2.*ProbUUbar;  // Probability of ssbar production 
    Prob_QQbar[3]=0.0;               // Probability of ccbar production
    Prob_QQbar[4]=0.0;               // Probability of bbbar production

    for ( G4int i=0 ; i<350 ; i++ ) { // Must be checked
      FS_LeftHadron[i] = 0;
      FS_RightHadron[i] = 0;
      FS_Weight[i] = 0.0; 
    }

    NumberOf_FS = 0;
}

// --------------------------------------------------------------

void G4VLongitudinalStringDecay::SetMinimalStringMass(const G4FragmentingString * const string)  
{
        //MaxMass = -350.0*GeV;
	G4double EstimatedMass=MaxMass;

        G4ParticleDefinition* LeftParton  = string->GetLeftParton();
        G4ParticleDefinition* RightParton = string->GetRightParton();
        if( LeftParton->GetParticleSubType() == RightParton->GetParticleSubType() ) { // q qbar, qq qqbar
          if( LeftParton->GetPDGEncoding() * RightParton->GetPDGEncoding() > 0 ) {
            // Not allowed combination of the partons
            throw G4HadronicException(__FILE__, __LINE__,
				      "G4VLongitudinalStringDecay::SetMinimalStringMass: Illegal quark content as input");
          }
        }
        if( LeftParton->GetParticleSubType() != RightParton->GetParticleSubType() ) { // q qq, qbar qqbar
          if( LeftParton->GetPDGEncoding() * RightParton->GetPDGEncoding() < 0 ) {
            // Not allowed combination of the partons
            throw G4HadronicException(__FILE__, __LINE__,
				      "G4VLongitudinalStringDecay::SetMinimalStringMass: Illegal quark content as input");
          }
        }
	
        G4int Qleft =std::abs(string->GetLeftParton()->GetPDGEncoding());
        G4int Qright=std::abs(string->GetRightParton()->GetPDGEncoding());

        if ((Qleft < 6) && (Qright < 6)) {   // Q-Qbar string
          EstimatedMass=minMassQQbarStr[Qleft-1][Qright-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if ((Qleft < 6) && (Qright > 1000)) {   // Q - DiQ string
          G4int q1=Qright/1000;
          G4int q2=(Qright/100)%10;
          EstimatedMass=minMassQDiQStr[Qleft-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;                    // It can be negative!
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if ((Qleft > 1000) && (Qright < 6)) {   // DiQ - Q string   6 6 6
          G4int q1=Qleft/1000;
          G4int q2=(Qleft/100)%10;
          EstimatedMass=minMassQDiQStr[Qright-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;                    // It can be negative!
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        // DiQuark - Anti DiQuark string -----------------

	G4double StringM=string->Get4Momentum().mag();

        #ifdef debug_LUNDfragmentation
        // G4cout<<"MinStringMass// Input String mass "<<string->Get4Momentum().mag()<<" Qleft "<<Qleft<<G4endl;
        #endif

        G4int q1= Qleft/1000    ;
        G4int q2=(Qleft/100)%10 ;

        G4int q3= Qright/1000   ;
        G4int q4=(Qright/100)%10;

        // -------------- 2 baryon production or 2 mesons production --------

        G4double EstimatedMass1 = minMassQDiQStr[q1-1][q2-1][0];
        G4double EstimatedMass2 = minMassQDiQStr[q3-1][q4-1][0];
        // Mass is negative if there is no corresponding particle.

        if ( (EstimatedMass1 > 0.) && (EstimatedMass2 > 0.)) {
           EstimatedMass = EstimatedMass1 + EstimatedMass2;
           if ( StringM > EstimatedMass ) {                     // 2 baryon production is possible.
              MinimalStringMass=EstimatedMass1 + EstimatedMass2;
              SetMinimalStringMass2(EstimatedMass);
              return;
          }
        }

        if ( (EstimatedMass1 < 0.) && (EstimatedMass2 > 0.)) {
           EstimatedMass = MaxMass;
           MinimalStringMass=EstimatedMass;
           SetMinimalStringMass2(EstimatedMass);
           return;
        }

        if ( (EstimatedMass1 > 0.) && (EstimatedMass2 < 0.)) {
           EstimatedMass = EstimatedMass1;
           MinimalStringMass=EstimatedMass;
           SetMinimalStringMass2(EstimatedMass);
           return;
        }

        //      if ( EstimatedMass >= StringM ) {
        // ------------- Re-orangement ---------------
        EstimatedMass=std::min(minMassQQbarStr[q1-1][q3-1] + minMassQQbarStr[q2-1][q4-1],
                               minMassQQbarStr[q1-1][q4-1] + minMassQQbarStr[q2-1][q3-1]);

        // In principle, re-arrangement and 2 baryon production can compete.
        // More physics consideration is needed.

        MinimalStringMass=EstimatedMass;
        SetMinimalStringMass2(EstimatedMass);

        return;
}

//--------------------------------------------------------------------------------------

void G4VLongitudinalStringDecay::SetMinimalStringMass2(const G4double aValue)
{
	MinimalStringMass2=aValue * aValue;
}

