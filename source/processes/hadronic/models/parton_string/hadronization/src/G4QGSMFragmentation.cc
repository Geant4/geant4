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

#include "G4Pow.hh"

//#define debug_QGSMfragmentation 

// Class G4QGSMFragmentation 
//****************************************************************************************
 
G4QGSMFragmentation::G4QGSMFragmentation() :
arho(0.5), aphi(0.), an(-0.5), ala(-0.75), aksi(-1.), alft(0.5)
{
    SetStrangenessSuppression((1.0 - 0.16)/2.);
    SetDiquarkSuppression(0.299);
    SetDiquarkBreakProbability(0.7);

    //... pspin_meson is probability to create pseudo-scalar meson 
    pspin_meson = 0.25;  SetVectorMesonProbability(pspin_meson);

    //... pspin_barion is probability to create 1/2 barion 
    pspin_barion = 0.5;  SetSpinThreeHalfBarionProbability(pspin_barion);

    //... vectorMesonMix[] is quark mixing parameters for vector mesons (Variable spin = 3)
    vectorMesonMix[0] = 0.;     //0.5 
    vectorMesonMix[1] = 0.375;  //0.5;
    vectorMesonMix[2] = 0.0;  
    vectorMesonMix[3] = 0.375;  //0.5; 
    vectorMesonMix[4] = 1.0;
    vectorMesonMix[5] = 1.0; 

    SetVectorMesonMixings(vectorMesonMix);
    SetMinMasses();
}

G4QGSMFragmentation::~G4QGSMFragmentation()
{}

//----------------------------------------------------------------------------------------------------------

G4KineticTrackVector* G4QGSMFragmentation::FragmentString(const G4ExcitedString& theString)
{
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
  G4KineticTrackVector * LeftVector=LightFragmentationTest(&theString);

  #ifdef debug_QGSMfragmentation
  if ( LeftVector != 0 ) G4cout<<"Non fragmentable - the string is converted to one hadron "<<G4endl;
  #endif

  if ( LeftVector != 0 ) return LeftVector;

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
                G4cout<<"Split remaining string into 2 final hadrons."<<G4endl;
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

G4double G4QGSMFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  
                                            G4ParticleDefinition* pHadron, G4double , G4double )
{    
  #ifdef debug_QGSMfragmentation
  G4cout<<"GetLightConeZ zmin zmax Parton pHadron "<<zmin<<" "<<zmax<<" "<<PartonEncoding<<" "<<pHadron->GetParticleName()<<G4endl;
  #endif

  G4double z;    
  G4double d1, d2, yf;
  G4double invD1(0.),invD2(0.), r1(0.),r2(0.),r12(0.);

  G4int absCode = std::abs( PartonEncoding );
  G4int absHadronCode=std::abs(pHadron->GetPDGEncoding());

  G4int q1, q2, q3;
  q1 = absHadronCode/1000; q2 = (absHadronCode % 1000)/100; q3 = (absHadronCode % 100)/10;

  G4bool StrangeHadron = (q1 == 3) || (q2 == 3) || (q3 == 3);

  if (absCode < 10)
  {                                              // A quark fragmentation ----------------------------
    if (absCode == 1 || absCode == 2) 
    {
      if (absHadronCode < 1000) 
      {                        // Meson  produced
        if ( !StrangeHadron )  {d1=2.0;        d2 = -arho + alft;}  // d1=2.0;
        else                   {d1=1.0;        d2 = -aphi + alft;}  // d1=1.0;
      } 
      else                 
      {                        // Baryon produced
        if ( !StrangeHadron )  {d1=0.0;        d2 =      arho - 2.0*an        + alft;}
        else                   {d1=0.0;        d2 =  2.0*arho - 2.0*an - aphi + alft;} 
      }
    } 
    else if (absCode == 3)  
    {
      if (absHadronCode < 1000) {d1=1.0 - aphi; d2 =  -arho          + alft;}  // Meson  produced s->K + u/d   d1=1.0
      else                      {d1=1.0 - aphi; d2 =   arho - 2.0*an + alft;}  // Baryon produced 

    } else throw G4HadronicException(__FILE__, __LINE__, "Unknown PDGencoding in G4QGSMFragmentation::G4LightConeZ");

    #ifdef debug_QGSMfragmentation
    G4cout<<"d1 d2 "<<d1<<" "<<d2<<G4endl;
    #endif

    d1+=1.0; d2+=1.0;

    invD1=1./d1; invD2=1./d2;

    const G4int maxNumberOfLoops = 10000;
    G4int loopCounter = 0;
    do
    {
      r1=G4Pow::GetInstance()->powA(G4UniformRand(),invD1);
      r2=G4Pow::GetInstance()->powA(G4UniformRand(),invD2);
      r12=r1+r2;
      z=r1/r12;
    } while ( ( (r12 > 1.0) || !((zmin <= z)&&(z <= zmax))) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
    }

    return z;
  }
  else
  {                                              // A di-quark fragmentation -------------------------
    if (absCode == 1103 || absCode == 2101 || 
        absCode == 2203 || absCode == 2103)
    {
      if(absHadronCode < 1000)                                                // Meson production
      {
         if ( !StrangeHadron )  {d1=1.0; d2=     arho - 2.0*an        + alft;} 
         else                   {d1=0.0; d2 = 2.*arho - 2.0*an - aphi + alft;} // d1=1.0;
      } 
      else                                                                 // Baryon production
      {
        if ( !StrangeHadron )  {d1=2.0*(arho - an); d2= -arho         + alft;}
        else                   {d1=2.0*(arho - an); d2 =-aphi         + alft;} 
      }

      #ifdef debug_QGSMfragmentation
      G4cout<<"d1 d2 "<<d1<<" "<<d2<<G4endl;
      #endif

      d1+=1.0; d2+=1.0;
      invD1=1./d1; invD2=1./d2;

      const G4int maxNumberOfLoops = 10000;
      G4int loopCounter = 0;
      do
      {
        r1=G4Pow::GetInstance()->powA(G4UniformRand(),invD1);
        r2=G4Pow::GetInstance()->powA(G4UniformRand(),invD2);
        r12=r1+r2;
        z=r1/r12;
      } while ( ( (r12 > 1.0) || !((zmin <= z)&&(z <= zmax))) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */ 
      if ( loopCounter >= maxNumberOfLoops ) {
        z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
      }

      return z;
    }
    else if (absCode == 3101 || absCode == 3103 ||           // For strange d-quarks
             absCode == 3201 || absCode == 3203)
    {
      // For future improvements
      // if (absHadronCode < 1000) {d1=1.0;}  // Meson production
      // else                      {d1=2.0;}  // Baryon production
      d2 =  (alft - (2.*ala - arho));
    }
    else
    {
      // if (absHadronCode < 1000) {d1=1.0;}      // Meson production
      //      else                     {d1=2.0;}  // Baryon production
      d2 =  (alft - (2.*aksi - arho));
    }

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do  
    {
      z = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      yf = G4Pow::GetInstance()->powA(d1, d2);
    } 
    while ( (G4UniformRand() > yf) && ++loopCounter < maxNumberOfLoops );  // Loop checking, 07.08.2015, A.Ribon //   
    /* For future improvements
    d1+=1.0; d2+=1.0;
    invD1=1./d1; invD2=1./d2;
    const G4int maxNumberOfLoops = 10000;
    G4int loopCounter = 0;
    do
    {
      r1=G4Pow::GetInstance()->powA(G4UniformRand(),invD1);
      r2=G4Pow::GetInstance()->powA(G4UniformRand(),invD2);
      r12=r1+r2;
      z=r1/r12;
    } while( ( (r12 > 1.0) || !((zmin <= z)&&(z <= zmax))) && ++loopCounter < maxNumberOfLoops );
    */
    if ( loopCounter >= maxNumberOfLoops ) {
       z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
    }

    return z;

  }

  return z;
}
//-----------------------------------------------------------------------------------------

G4LorentzVector * G4QGSMFragmentation::SplitEandP(G4ParticleDefinition * pHadron,
                                                  G4FragmentingString * string,  
                                                  G4FragmentingString * NewString)
{
       G4double HadronMass = pHadron->GetPDGMass();

       SetMinimalStringMass(NewString);

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

       //...  sample Pt of the hadron
       G4int attempt=0;
       do
       {
         attempt++; if (attempt > StringLoopInterrupt) return 0;

         HadronPt =SampleQuarkPt()  + string->DecayPt();	
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
    G4ParticleDefinition * LeftHadron, * RightHadron;
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
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

	  G4int IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? -1 : +1; // if we have a quark, 
                                                                          // we need antiquark or diquark
	  pDefPair QuarkPair = CreatePartonPair(IsParticle);
	  quark = QuarkPair.second;

	  LeftHadron= hadronizer->BuildLowSpin(QuarkPair.first, string->GetLeftParton());
	} else {
	  //... there is a Diquark on cluster ends
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
      	  LeftHadron=hadronizer->BuildLowSpin(QuarkPair.first, string->GetLeftParton());
        }

        RightHadron = hadronizer->BuildLowSpin(string->GetRightParton(), quark);

    } while ( ( ResidualMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass() )
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

G4bool G4QGSMFragmentation::IsFragmentable(const G4FragmentingString * const string)
{
	return sqr( FragmentationMass(string) + MassCut ) < string->Mass2();
}

//----------------------------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::StopFragmenting(const G4FragmentingString * const string)
{
	SetMinimalStringMass(string);
	if (string->FourQuarkString())
	{
	  return G4UniformRand() < G4Exp(-0.0005*(string->Mass() - MinimalStringMass));
	} else { 
	  G4bool Result = G4UniformRand() < 
			  G4Exp(-0.66e-6*(string->Mass()*string->Mass() - MinimalStringMass*MinimalStringMass));
          // G4bool Result = string->Mass() < MinimalStringMass + 150.*MeV*G4UniformRand();  // a'la LUND

          #ifdef debug_QGSMfragmentation 
          G4cout<<"StopFragmenting MinimalStringMass string->Mass() "<<MinimalStringMass<<" "<<string->Mass()<<G4endl; 
          G4cout<<"StopFragmenting - Yes/No "<<Result<<G4endl;
          #endif	
	  return Result;
	}
}


//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::Sample4Momentum(G4LorentzVector* Mom    , G4double Mass    , 
                                          G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass) 
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
      SetStrangenessSuppression((1.0 - 0.07)/2.);
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

      return had;

   } else {  //... Diquark does not break

      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1;  // if we have a diquark, we need quark)
      G4double StrSup=GetStrangeSuppress();  // for changing s-sbar production
      SetStrangenessSuppression((1.0 - 0.07)/2.);
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      SetStrangenessSuppression(StrSup);

      created = QuarkPair.second;

      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
      return had;
   }
}

//-----------------------------------------------------------------------------

G4KineticTrack * G4QGSMFragmentation::Splitup( G4FragmentingString *string, 
			                       G4FragmentingString *&newString)
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

//---------------------------------------------------------------
void G4QGSMFragmentation::SetMinMasses()
{
    // ------ For estimation of a minimal string mass ---------------
    Mass_of_light_quark    =140.*MeV;
    Mass_of_heavy_quark    =500.*MeV;
    Mass_of_string_junction=720.*MeV;

    G4double minMQQbarStr[3][3] ={ {350.0, 350.0, 710.0},  //DDbar, DUbar, DSbar in MeV
                                   {350.0, 350.0, 710.0},  //UDbar, UUbar, USbar in Mev
                                   {710.0, 710.0,1070.0 }};//SDbar, SUbar, SSbar in MeV
    for (G4int i=0; i<3; i++){ for (G4int j=0; j<3; j++){minMassQQbarStr[i][j]=minMQQbarStr[i][j];};}; 

 
    G4double minMQDiQStr[3][3][3] = {{{1160., 1160., 1340.}, {1160., 1160., 1340.}, {1340., 1340., 1540.},}, //d-dd, d-du, d-ds, d-ud, d-uu, d-us, d-sd, d-su, d-ss 
                                    {{1160., 1160., 1340.}, {1160., 1160., 1340.}, {1340., 1340., 1540.},}, //u-dd, u-du, u-ds, u-ud, u-uu, u-us, u-sd, u-su, u-ss
                                    {{1520., 1520., 1690.}, {1520., 1520., 1690.}, {1690., 1690., 1890. }}};//s-dd, s-du, s-ds, s-ud, s-uu, s-us, s-sd, s-su, s-ss
    for (G4int i=0; i<3; i++){ for (G4int j=0; j<3; j++){ for (G4int k=0; k<3; k++){minMassQDiQStr[i][j][k]=minMQDiQStr[i][j][k];};};};                      

    // ------ An estimated minimal string mass ----------------------
    MinimalStringMass  = 0.;
    MinimalStringMass2 = 0.;
}

//--------------------------------------------------------------------------------------
void G4QGSMFragmentation::SetMinimalStringMass(const G4FragmentingString  * const string)  
{
	G4double EstimatedMass=0.;

        G4int Qleft =std::abs(string->GetLeftParton()->GetPDGEncoding());
        G4int Qright=std::abs(string->GetRightParton()->GetPDGEncoding());

        if ((Qleft < 4) && (Qright < 4)) {   // Q-Qbar string
          EstimatedMass=minMassQQbarStr[Qleft-1][Qright-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if ((Qleft < 4) && (Qright > 1000)) {   // Q - DiQ string
          G4int q1=Qright/1000;
          G4int q2=(Qright/100)%10;
          EstimatedMass=minMassQDiQStr[Qleft-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        if ((Qleft > 1000) && (Qright < 4)) {   // DiQ - Q string
          G4int q1=Qleft/1000;
          G4int q2=(Qleft/100)%10;
          EstimatedMass=minMassQDiQStr[Qright-1][q1-1][q2-1];
          MinimalStringMass=EstimatedMass;
          SetMinimalStringMass2(EstimatedMass);
          return;
        }

        // DiQuark - Anti DiQuark string -----------------
	G4int Number_of_quarks=0;
        G4int Number_of_squarks=0;
        
	G4double StringM=string->Get4Momentum().mag();

        #ifdef debug_QGSMfragmentation
        // G4cout<<"MinStringMass// Input String mass "<<string->Get4Momentum().mag()<<" Qleft "<<Qleft<<G4endl;
        #endif

	if ( Qleft > 1000)
	{
		Number_of_quarks+=2;
		G4int q1=Qleft/1000;
		if ( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if ( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

		G4int q2=(Qleft/100)%10;
		if ( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if ( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
	}

        #ifdef debug_QGSMfragmentation
        // G4cout<<"Min mass with Qleft "<<Qleft<<" "<<EstimatedMass<<G4endl;
        #endif

	if ( Qright > 1000)
	{
		Number_of_quarks+=2;
		G4int q1=Qright/1000;
		if ( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if ( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

		G4int q2=(Qright/100)%10;
		if ( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
		if ( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
		//EstimatedMass +=Mass_of_string_junction;
	}

        #ifdef debug_QGSMfragmentation
        // G4cout<<"Min mass with Qleft and Qright "<<Qright<<" "<<EstimatedMass<<G4endl;
        // G4cout<<"Number_of_quarks "<<Number_of_quarks<<" Number_of_squarks "<<Number_of_squarks<<G4endl;
        #endif

	if (Number_of_quarks==4)
	{
		if (StringM > 1880.) {                                          // 2*Mn = 1880
          	  if (Number_of_squarks==0)      {EstimatedMass += 1320.*MeV;}  //560+1320=1880=2*Mn
    	          else if (Number_of_squarks==1) {EstimatedMass += 1150.*MeV;}  //920+1150=2070=M(Lam+N)
          	  else if (Number_of_squarks==2) {EstimatedMass +=  960.*MeV;}  //1280+960=2240= 2*M Lam
          	  else if (Number_of_squarks==3) {EstimatedMass +=  800.*MeV;}  //1640+800=2440=Mxi+Mlam
          	  else if (Number_of_squarks==4) {EstimatedMass +=  640.*MeV;}  //2000+640=2640=2*Mxi
                  else {}
                }
		else
		{
          	  if (Number_of_squarks < 3)     {EstimatedMass -= 200.*MeV;}
          	  else if (Number_of_squarks==3) {EstimatedMass -=  50.*MeV;}
          	  else if (Number_of_squarks==4) {EstimatedMass -=  40.*MeV;}
                  else {}
		}
	}

        #ifdef debug_QGSMfragmentation
        // G4cout<<"EstimatedMass -------------------- "<<EstimatedMass <<G4endl;
        #endif

	MinimalStringMass=EstimatedMass;
	SetMinimalStringMass2(EstimatedMass);
}

//--------------------------------------------------------------------------------------
void G4QGSMFragmentation::SetMinimalStringMass2(const G4double aValue)
{
	MinimalStringMass2=aValue * aValue;
}

