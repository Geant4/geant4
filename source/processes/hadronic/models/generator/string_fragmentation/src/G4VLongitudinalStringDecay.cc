// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VLongitudinalStringDecay.cc,v 1.6 1999/12/15 14:52:48 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//  Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
// -----------------------------------------------------------------------------
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4VLongitudinalStringDecay.hh"

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
   
   StrangeSuppress  = 0.44;
   DiquarkSuppress  = 0.1;
   DiquarkBreakProb = 0.1;
   
   //... pspin_meson is probability to create vector meson 
   pspin_meson = 0.5;

   //... pspin_barion is probability to create 3/2 barion 
   pspin_barion = 0.5;

   //... pmix_meson0[] is quark mixing parameters for scalar mesons (Variable spin=1)
   pmix_meson1[0] = 0.5;
   pmix_meson1[1] = 0.0;
   pmix_meson1[2] = 0.5;
   pmix_meson1[3] = 0.0;
   pmix_meson1[4] = 1.0;
   pmix_meson1[5] = 1.0; 

   //... pmix_meson1[] is quark mixing parameters for vector mesons (Variable spin = 3)
   pmix_meson0[0] = 0.5; 
   pmix_meson0[1] = 0.25; 
   pmix_meson0[2] = 0.5; 
   pmix_meson0[3] = 0.25; 
   pmix_meson0[4] = 1.0; 
   pmix_meson0[5] = 0.5; 

// Parameters may be changed until the first fragmentation starts
   PastInitPhase=false;
   
}
   

G4VLongitudinalStringDecay::~G4VLongitudinalStringDecay()
   {
   }

//********************************************************************************
// Operators

//const  & G4VLongitudinalStringDecay::operator=(const G4VLongitudinalStringDecay &right)
//    {
//    }

int G4VLongitudinalStringDecay::operator==(const G4VLongitudinalStringDecay &right) const
    {
	G4Exception("G4VLongitudinalStringDecay::operator== forbidden");
	return false;
//    return !memcmp(this, &right, sizeof(G4VLongitudinalStringDecay));
    }

int G4VLongitudinalStringDecay::operator!=(const G4VLongitudinalStringDecay &right) const
    {
	G4Exception("G4VLongitudinalStringDecay::operator!= forbidden");
	return true;
//    return memcmp(this, &right, sizeof(G4VLongitudinalStringDecay));
    }


//********************************************************************************

G4int G4VLongitudinalStringDecay::SampleQuarkFlavor(void)
   {
   return (1 + (int)(G4UniformRand()/StrangeSuppress));
   }

//********************************************************************************

void G4VLongitudinalStringDecay::SampleQuarkPt(G4double* thePx, G4double* thePy)
   {
   G4double Pt = -log(G4UniformRand());
   Pt = sqrt(Pt);
   Pt *= SigmaQT;
   G4double phi = 2.*pi*G4UniformRand();
   *thePx = Pt * cos(phi);
   *thePy = Pt * sin(phi);
   }


//********************************************************************************

void G4VLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector* Hadrons)
   {
   // `yo-yo` formation time
   const G4double kappa = 1.0 * GeV/fermi;
   for(G4int c1 = 0; c1 < Hadrons->length(); c1++)
      {
      G4double SumPz = 0; 
      G4double SumE  = 0;
      for(G4int c2 = 0; c2 < c1; c2++)
         {
         SumPz += Hadrons->at(c2)->Get4Momentum().pz();
         SumE  += Hadrons->at(c2)->Get4Momentum().e();   
         } 
      G4double HadronE  = Hadrons->at(c1)->Get4Momentum().e();
      G4double HadronPz = Hadrons->at(c1)->Get4Momentum().pz();
      Hadrons->at(c1)->SetFormationTime((theInitialStringMass - 2.*SumPz + HadronE - HadronPz)/(2.*kappa));
      G4ThreeVector aPosition(0, 0,     (theInitialStringMass - 2.*SumE  - HadronE + HadronPz)/(2.*kappa));
      Hadrons->at(c1)->SetPosition(aPosition);
      } 
   }

//********************************************************************************
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
//********************************************************************************

G4ParticleDefinition* G4VLongitudinalStringDecay::CreateHadron(G4int id1, G4int id2, G4bool theGivenSpin, G4int theSpin) 
   {

   G4int PDGEncoding;
   if (abs(id1) < abs(id2))
      {
      int xchg = id1; 
      id1 = id2;  
      id2 = xchg;
      }       
   G4int ifl1 = abs(id1);
   //...  Construct meson with account flavor mixing     
   if (abs(id1) < 1000 && abs(id2) <  1000) 
      {
      G4int spin = theSpin;
      if(!theGivenSpin) 
         spin = 2*(G4int)(pspin_meson + G4UniformRand()) + 1;
      if (id1 + id2 == 0) 
         {    
         G4double rmix = G4UniformRand();
         G4int    imix = 2*ifl1 - 1;
         if(spin == 1)
            {
            PDGEncoding = 110*(1 + (G4int)(rmix + pmix_meson0[imix - 1])
               + (G4int)(rmix + pmix_meson0[imix])) +  spin;
            }
         else
            {
            PDGEncoding = 110*(1 + (G4int)(rmix + pmix_meson1[imix - 1])
               + (G4int)(rmix + pmix_meson1[imix])) +  spin;
            }
         }
      else
         {
         G4int kfs = abs((id2 < 0)? id2 : id1);
         PDGEncoding = 100 * ifl1 + 10 * abs(id2) +  spin;  
         if((!(ifl1&1) && kfs == ifl1) || ((ifl1&1) && kfs != ifl1)) 
            PDGEncoding = - PDGEncoding;
         }
      }
   else
      {
      ifl1 /= 1000;
      G4int ifl2 = (abs(id1) - ifl1 * 1000)/100;
      G4int kflds = abs(id1)%10; 
      G4int ifl3 = id2;
      if (id1 < 0)
         {
         ifl1 = - ifl1;
         ifl2 = - ifl2;
         }
      //... Construct barion, distinguish Lambda and Sigma barions.
      G4int kfla = abs(ifl1);
      G4int kflb = abs(ifl2);
      G4int kflc = abs(ifl3);

      G4int kfld = G4std::max(kfla,kflb);
            kfld = G4std::max(kfld,kflc);
      G4int kflf = G4std::min(kfla,kflb);
            kflf = G4std::min(kflf,kflc);

      G4int kfle = kfla + kflb + kflc - kfld - kflf;
      //... barion with content uuu or ddd or sss has always spin = 4
      G4int spin = (kfla == kflb && kflb == kflc)? 4: theSpin;   
      if(!theGivenSpin)
	 spin = (pspin_barion > G4UniformRand() || (kfla == kflb &&  kfla == kflc))? 4 : 2;

      G4int kfll = 0;
      if(spin == 2 && kfld > kfle && kfle > kflf)
         {
         if(kflds == 1 && kfla == kfld)
             kfll = 1;
         if(kflds == 1 && kfla != kfld)
             kfll = (G4int)(0.25 + G4UniformRand());
         if(kflds == 3 && kfla != kfld)
             kfll = (G4int)(0.75 + G4UniformRand());
         }
      if (kfll == 1)
	 PDGEncoding = 1000 * kfld + 100 * kflf + 10 * kfle + spin;
      else    
	 PDGEncoding = 1000 * kfld + 100 * kfle + 10 * kflf + spin;
      if (id1 < 0)
	 PDGEncoding = -PDGEncoding;
      }
   return FindParticle(PDGEncoding); 
   }
    
//*******************************************************************************************************
//*******************************************************************************

void G4VLongitudinalStringDecay::SimpleString::update()
{
	MassSquare = left->w()*right->w() - sqr(left->Px() + right->Px()) - sqr(left->Py() + right->Py());
}
G4ParticleDefinition * G4VLongitudinalStringDecay::SimpleString::Splitup(G4int & NewDecayEncoding)  
{
      //...the main iteration loop for the fragmentation is started here:   
      //... stable and non-stable hadrons can be produced ... update theGivenSpin flag
      //...theGivenSpin = false;
      //...Initialize di_quark_break flag



       //... random choice of string end to use for creating the hadron (decay)   
       Side = (G4UniformRand() < 0.5)? 1: -1;
       if (Side < 0)
       {
          decay  = Right(); 
          stable = Left();
       } else
       {
          decay  = Left();
	  stable = Right();
       }
	  
       //... if string end is  a quark  
       G4int signdecayEncoding = decay->Encoding()/abs(decay->Encoding());
       G4int QuarkEncoding;
       if (abs(decay->Encoding()) < 1000) 
          {    
          //... it is quark ...  Generate q,qbar or qq,qqbar pair    
          if (G4UniformRand() >= theStringDecay->GetDiquarkSuppress())
             {
             //... q,qbar pair is choosen, sample quark flavor ifln 
             QuarkEncoding    = -signdecayEncoding*theStringDecay->SampleQuarkFlavor();
             NewDecayEncoding = -QuarkEncoding;
             }
          else
             {
             //... sample quarks and Construct Diquark ifln 
             G4int i1  = theStringDecay->SampleQuarkFlavor();
             G4int i2  = theStringDecay->SampleQuarkFlavor();
             G4int i10 = G4std::max(i1,i2);
             G4int i20 = G4std::min(i1,i2);
             G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
             QuarkEncoding = signdecayEncoding * (i10 * 1000 + i20 * 100 + spin);
             NewDecayEncoding = -QuarkEncoding;
             }    
          } 
       else   
          { 
	  //... it is Diquark ...  
	  //... can Diquark break or not? 
          if (G4UniformRand() < theStringDecay->GetDiquarkBreakProb())
             {
             //... Diquark break

	     G4int stableQuarkEncoding = decay->Encoding()/1000;
	     G4int decayQuarkEncoding = (decay->Encoding()/100)%10;
	     if (G4UniformRand() < 0.5)
                {
                G4int Swap = stableQuarkEncoding;
                stableQuarkEncoding = decayQuarkEncoding;
                decayQuarkEncoding = Swap;
                }
	     Decay()->setEncoding(decayQuarkEncoding);
	     	
             QuarkEncoding = -decayQuarkEncoding/abs(decayQuarkEncoding)*theStringDecay->SampleQuarkFlavor();
             //... Build new Diquark
             G4int i10  = G4std::max(abs(QuarkEncoding), abs(stableQuarkEncoding));
             G4int i20  = G4std::min(abs(QuarkEncoding), abs(stableQuarkEncoding));
             G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
             NewDecayEncoding = -QuarkEncoding/abs(QuarkEncoding)*(i10 * 1000 + i20 * 100 + spin);
             }
          else
             {
             //... Diquark does not break 
             QuarkEncoding = signdecayEncoding*theStringDecay->SampleQuarkFlavor();
             NewDecayEncoding = -QuarkEncoding;
             }
          }
       //...  Construct produced hadron: encoding, mass and transverse momentum 

       G4ParticleDefinition * pHadron = theStringDecay->CreateHadron(decay->Encoding(), QuarkEncoding, false,  1);
       return pHadron;
}

G4bool G4VLongitudinalStringDecay::SimpleString::SplitLast(G4KineticTrackVector * LeftVector,
							   G4KineticTrackVector * RightVector)
{
    //... perform last cluster decay
    G4double ResidualStringPx  = Left()->Px() + Right()->Px();
    G4double ResidualStringPy  = Left()->Py() + Right()->Py();
    G4double ResidualMass     = sqrt(MassSquared());
    G4double ResidualStringPz = (Left()->w() - Right()->w())/2.;
    G4double ResidualStringE  = (Left()->w() + Right()->w())/2.;
    G4ThreeVector ClusterVel(ResidualStringPx/ResidualStringE, ResidualStringPy/ResidualStringE, ResidualStringPz/ResidualStringE);
    G4ParticleDefinition* pLastHadron;
    G4double ClusterMassCut = theStringDecay->GetClusterMass();
    G4int cClusterInterrupt = 0;
    G4ParticleDefinition * LeftHadron, * RightHadron;
    do
       {
       if (cClusterInterrupt++ >= theStringDecay->GetClusterLoopInterrupt())
       {
          return false;
       }
       //... if any ifl is a Diquark 
       G4int signLeftEncoding = Left()->Encoding()/abs(Left()->Encoding());
       G4int QuarkEncoding;
       if (!(abs(Left()->Encoding()) > 1000  || abs(Right()->Encoding()) > 1000)) 
          {
	  //... there are quarks on cluster ends
	  //... randomly choose q,qbar pair or qq,qqbar pair
          if (G4UniformRand() > theStringDecay->GetDiquarkSuppress()) 
             {
             //... q,qbar pair is choosen, sample quark flavor ifln 
             QuarkEncoding = -signLeftEncoding*theStringDecay->SampleQuarkFlavor();
             }    
          else 
             {   
	     //... sample quarks and Construct Diquark ifln 
	     G4int i1 = theStringDecay->SampleQuarkFlavor();
	     G4int i2 = theStringDecay->SampleQuarkFlavor();
	     G4int i10 = G4std::max(i1,i2);
	     G4int i20 = G4std::min(i1,i2);
	     G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1: 3;
	     QuarkEncoding = signLeftEncoding*(i10 * 1000 + i20 * 100 + spin);
             }
          }
       else               
          {
	  //... there is a Diquark on cluster ends
	  //...  randomly choose q,qbar pair  
	  QuarkEncoding = ((abs(Left()->Encoding())  < 1000)? -signLeftEncoding : signLeftEncoding)*theStringDecay->SampleQuarkFlavor();
          }
       //...  encodings and masses of hadrons  
       //...  theGivenSpin = false; //means any hadron (with given quark content) can be built
       //...  theSpin = 1;
       LeftHadron =  theStringDecay->CreateHadron(Left()->Encoding(),  QuarkEncoding, false,  1);
       RightHadron = theStringDecay->CreateHadron(Right()->Encoding(), -QuarkEncoding, false,  1);

       //... repeat procedure, if mass of cluster is too low to produce hadrons
       //... ClusterMassCut = 0.15*GeV model parameter
       if (QuarkEncoding < 3)
               ClusterMassCut = 0.;
       } 
    while (ResidualMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()  + ClusterMassCut);

    //... compute hadron momenta and energies   
    G4LorentzVector  LeftMom, RightMom;
    G4ThreeVector    Pos;
    theStringDecay->Sample4Momentum(&LeftMom, LeftHadron->GetPDGMass(), &RightMom, RightHadron->GetPDGMass(), ResidualMass);
    LeftMom.boost(ClusterVel);
    RightMom.boost(ClusterVel);
    LeftVector->insert(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
    RightVector->insert(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

    return true;

}

//*****************************************************************************************************
//----------------------------------------------------------------------------------------------------------
G4KineticTrackVector* G4VLongitudinalStringDecay::FragmentString(const G4ExcitedString& theString)
{
//    Can no longer modify Parameters for Fragmentation.
	PastInitPhase=true;
	
// 	check if string has enough mass to fragment...
	G4KineticTrackVector * LeftVector=LightFragmentationTest(theString);
	if ( LeftVector != 0 ) return LeftVector;
	
	LeftVector = new G4KineticTrackVector;
	G4KineticTrackVector * RightVector=new G4KineticTrackVector;
	
	G4bool success=false;
	
	for ( G4int attempts=0; attempts < StringLoopInterrupt; attempts++)
	{
		SimpleString currentString(theString, this);
		LeftVector->clearAndDestroy();
		RightVector->clearAndDestroy();

		while (! StopFragmenting(currentString) )
		{  // Split current string into hadron + new string
			G4int newDecayEncoding;  // used as output from SplitUp...
			G4ParticleDefinition * hadronDefinition=currentString.Splitup(newDecayEncoding);
			G4KineticTrack * Hadron=SplitEandP(hadronDefinition,&currentString, newDecayEncoding);
			if ( Hadron == 0 ) goto ENDLOOP;
			if ( currentString.GetDecayDirection() > 0 )
				LeftVector->insert(Hadron);
       			else
	  			RightVector->insert(Hadron);
			currentString.update();
			if ( ! IsFragmentable(currentString) ) goto ENDLOOP;
		} 
		// Split current string into 2 final Hadrons
		if ( ! currentString.SplitLast(LeftVector, RightVector) ) goto ENDLOOP;
		success=true;
		break;   // Success!!!
ENDLOOP:	;	
		
	}
	
	if ( ! success )
	{
		LeftVector->clearAndDestroy();
		RightVector->clearAndDestroy();
		return LeftVector;
	}
		
	// After the next procedure LeftVector will be contained List of Hadrons
	//	  in correct order.
	while(!RightVector->isEmpty())
	    LeftVector->insert(RightVector->removeLast());
	delete RightVector;

	CalculateHadronTimePosition(theString.Get4Momentum().mag(), LeftVector);

	// Make backward lorenz boost
	G4LorentzRotation toCms(-1*theString.Get4Momentum().boostVector());
	G4LorentzVector pLeftRest= toCms*theString.GetLeftParton()->Get4Momentum();
	G4LorentzVector pRightRest= toCms*theString.GetRightParton()->Get4Momentum();

	toCms.rotateZ(-1*pLeftRest.phi());
	toCms.rotateY(-1*pLeftRest.theta());

	pLeftRest=toCms*theString.GetLeftParton()->Get4Momentum();
	pRightRest=toCms*theString.GetRightParton()->Get4Momentum();
	G4LorentzRotation toObserverFrame(toCms.inverse());

	for(int C1 = 0; C1 < LeftVector->length(); C1++)
	{
	   G4KineticTrack* Hadron = LeftVector->at(C1);
	   G4LorentzVector Momentum = Hadron->Get4Momentum();
	   Momentum = toObserverFrame*Momentum;
	   Hadron->Set4Momentum(Momentum);
	   G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
	   Momentum = toObserverFrame*Coordinate;
	   Hadron->SetFormationTime(Momentum.e());
	   G4ThreeVector aPosition(Momentum.px(), Momentum.py(), Momentum.pz());
	   Hadron->SetPosition(theString.GetPosition()+aPosition);
	}
	return LeftVector;
		


}
G4KineticTrack * G4VLongitudinalStringDecay::SplitEandP(G4ParticleDefinition * pHadron,
	G4VLongitudinalStringDecay::SimpleString * string, G4int newDecayEncoding)
{
       G4double HadronMass = pHadron->GetPDGMass();

       // calculate and assign hadron transverse momentum component HadronPx andHadronPy
       G4double thePx, thePy;
       SampleQuarkPt(&thePx, &thePy);
       G4double HadronPx = string->Decay()->Px() + thePx;
       G4double HadronPy = string->Decay()->Py() + thePy;
       //...  sample z to define hadron longitudinal momentum and energy
       //... but first check the available phase space
       G4double DecayQuarkMass2  = sqr(FindParticle(string->Decay()->Encoding())->GetPDGMass());
       G4double HadronMass2T = HadronMass*HadronMass + HadronPx*HadronPx + HadronPy*HadronPy;
       if (DecayQuarkMass2 + HadronMass2T >= SmoothParam*string->Decay()->w()*string->Stable()->w()) 
          return 0;		// have to start all over!

       //... then compute allowed z region  z_min <= z <= z_max 
 
       G4double zMin = HadronMass2T/(string->Decay()->w()*string->Stable()->w());
       G4double zMax = 1. - DecayQuarkMass2/(string->Decay()->w()*string->Stable()->w());
       if (zMin >= zMax) return 0;		// have to start all over!
	
       
//####wasGF       G4double z = GetLightConeZ(zMin, zMax, QuarkEncoding, pHadron, HadronPx, HadronPy);      
       G4double z = GetLightConeZ(zMin, zMax, string->Decay()->Encoding(), pHadron, HadronPx, HadronPy);      
       
       //... now compute hadron longitudinal momentum and energy
       //longitudinal hadron momentum component HadronPz

       G4double HadronPz = (z * string->Decay()->w() - HadronMass2T/(z * string->Decay()->w()))*0.5;
       HadronPz *= string->GetDecayDirection();
       //total hadron energy HadronE
       G4double HadronE  = (z * string->Decay()->w() + HadronMass2T/(z * string->Decay()->w()))*0.5;

       //...update  (after hadron separation) string light-cone variables
       	
       string->Decay()->setEncoding(newDecayEncoding);
       string->Decay()->setpxpy(-thePx,-thePy);
       string->Left()->decreasew(HadronE + HadronPz);
       string->Right()->decreasew(HadronE - HadronPz);

       G4LorentzVector a4Momentum(HadronPx,HadronPy,HadronPz,HadronE);
       G4ThreeVector   Pos;
       G4KineticTrack * Hadron = new G4KineticTrack(pHadron, 0, Pos, a4Momentum);
       
       return Hadron;
}

G4double G4VLongitudinalStringDecay::MinFragmentationMass(SimpleString & theString,
					G4ParticleDefinition *& Hadron1,G4ParticleDefinition *& Hadron2)
{
        G4int LeftEncoding  = theString.Left()->Encoding();
        G4int RightEncoding = theString.Right()->Encoding();

        G4double mass;

        //      qq - (qq)bar
        G4bool FourQuarkString  = abs(LeftEncoding) > 1000 && abs(RightEncoding) > 1000;
        //      q - qbar
        G4bool QQbarString      = abs(LeftEncoding) < 1000 && abs(RightEncoding) < 1000;

        Hadron2=0;

        if (!FourQuarkString )
        {
           // spin 0 meson or spin 1/2 barion will be built

           G4int theSpinLight = QQbarString ? 1: 2;
           Hadron1 = CreateHadron(LeftEncoding, RightEncoding, true, theSpinLight);
           mass= (Hadron1)->GetPDGMass();
        } else
        {
           //... string is qq--qqbar type: Build two stable hadrons,
           //... but we need extra uubar or ddbar quark pair
           G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
           if (LeftEncoding < 0) iflc = -iflc;

           //... theSpin = 2; spin 1/2 baryons will be built
           Hadron1 = CreateHadron(LeftEncoding,  iflc, true, 2);
           Hadron2 = CreateHadron(RightEncoding,-iflc, true, 2);

           mass= (Hadron1)->GetPDGMass() + (Hadron2)->GetPDGMass();
        }
        return mass;
}
					
G4double G4VLongitudinalStringDecay::IsFragmentable(SimpleString & theString)
{
	G4ParticleDefinition * Hadron1=0, * Hadron2=0;

	return sqr(MinFragmentationMass(theString, Hadron1, Hadron2)+MassCut) < theString.MassSquared();
}

 
G4KineticTrackVector* G4VLongitudinalStringDecay::LightFragmentationTest(const G4ExcitedString& theString)
{
   // Check string decay threshold
	
	SimpleString aString(theString,this);
	
	G4KineticTrackVector * result=0;  // return 0 when string exceeds the mass cut -);
	
	G4ParticleDefinition * Hadron1=0, * Hadron2 =0;
	if ( sqr(MinFragmentationMass(aString, Hadron1, Hadron2)+MassCut) < theString.Get4Momentum().mag2())
	 return 0;
	
	result=new G4KineticTrackVector;
        
	if ( Hadron2 ==0 )
	{
	      	 // Substitute string by light hadron, Note that Energy is not conserved here!

	       G4ThreeVector Mom3 = theString.Get4Momentum().vect();
	       G4LorentzVector Mom(Mom3, sqrt(Mom3.mag2() + sqr(Hadron1->GetPDGMass())));
               result->insert(new G4KineticTrack(Hadron1, 0, theString.GetPosition(), Mom));
	} else 
	{
	   //... string was qq--qqbar type: Build two stable hadrons,
	       G4LorentzVector  Mom1, Mom2;
	       Sample4Momentum(&Mom1, Hadron1->GetPDGMass(), &Mom2, Hadron2->GetPDGMass(), theString.Get4Momentum().mag());
	       result->insert(new G4KineticTrack(Hadron1, 0, theString.GetPosition(), Mom1));
	       result->insert(new G4KineticTrack(Hadron2, 0, theString.GetPosition(), Mom2));
               G4ThreeVector Velocity = theString.Get4Momentum().boostVector();
               result->Boost(Velocity);          
	}
	return result;
	
}

G4bool G4VLongitudinalStringDecay::StopFragmenting(SimpleString& string)
{
	G4double SumMass = MassCut;
	
	//      qq - (qq)bar
	G4bool FourQuarkString	= abs(string.Left()->Encoding()) > 1000 && abs(string.Right()->Encoding()) > 1000;
	//      q - qbar 
	G4bool QQbarString	= abs(string.Left()->Encoding()) < 1000 && abs(string.Right()->Encoding()) < 1000;
 
	if ( ! FourQuarkString )
	{
	   //... string is q --qbar or q--qq type: Build a stable hadron
	   // spin 1 meson or spin 3/2 barion will be built
	   G4int theSpin = QQbarString ? 3: 4;	// 2*J+1 
	   SumMass+= CreateHadron(string.Left()->Encoding(), string.Right()->Encoding(), true, theSpin)->GetPDGMass();
	}else
	{
	   //... string is qq--qqbar type: Build two stable hadrons,
	      //   need extra uubar or ddbar quark pair
	   G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
	   if (string.Left()->Encoding() < 0) iflc = -iflc;

	   //... theSpin = 4; spin 3/2 baryons will be built
	   SumMass += CreateHadron(string.Left()->Encoding(),  iflc, true, 4)->GetPDGMass() +
	   	      CreateHadron(string.Right()->Encoding(),-iflc, true, 4)->GetPDGMass();
	}
	return string.MassSquared() <= sqr(SumMass); 
}


    
//*****************************************************************************************************

G4ParticleDefinition* G4VLongitudinalStringDecay::FindParticle(G4int Encoding) 
   {
   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);
   if (ptr == NULL) ptr = G4ParticleTable::GetParticleTable()->FindParticle(-Encoding);
   if (ptr == NULL)
       {
       G4cout << "Particle with encoding "<<Encoding<<" does not exist!!!"<<G4endl;
       G4Exception("Check your particle table");
       }
   return ptr;    
   }

//*****************************************************************************************************

G4bool G4VLongitudinalStringDecay::FragmentString(G4KineticTrackVector* aHadrons, const G4ExcitedString* theString)
   {
   G4KineticTrackVector* Output = FragmentString(*theString);
   if (Output)
       {
       while(!Output->isEmpty())
           aHadrons->insert(Output->removeLast());
       delete Output;
       return TRUE;
       }
   return FALSE;
   }


//*****************************************************************************************************


G4KineticTrackVector* G4VLongitudinalStringDecay::DecayResonans (G4KineticTrackVector* aHadrons)
    {
    G4KineticTrackVector* Out = new G4KineticTrackVector;
    for(G4int c1 = 0; c1 < aHadrons->length(); c1++)
       {
       G4KineticTrack* pKT = aHadrons->at(c1);
       if (abs(pKT->GetDefinition()->GetPDGEncoding())%10 > 2)
           {
      /*     G4KineticTrackVector* pKTV = pKT->Decay();
           if (pKTV)
              {
              while(!pKTV->isEmpty())
                Out->insert(pKTV->removeLast());
              pKTV->clearAndDestroy();
              delete pKTV;
              continue;
             }*/
           }
       Out->insert(new G4KineticTrack(*pKT));
       }
    return Out;
    }


//*****************************************************************************************************

void G4VLongitudinalStringDecay::Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass) 
    {
    G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) - sqr(2.*Mass*AntiMass);
    G4double Pabs = (r_val > 0.)? sqrt(r_val)/(2.*InitialMass) : 0;

    //... sample unit vector       
    G4double pz = 1. - 2.*G4UniformRand();  
    G4double st     = sqrt(1. - pz * pz)*Pabs;
    G4double phi    = 2.*pi*G4UniformRand();
    G4double px = st*cos(phi);
    G4double py = st*sin(phi);
    pz *= Pabs;
    
    Mom->setPx(px); Mom->setPy(py); Mom->setPz(pz);
    Mom->setE(sqrt(Pabs*Pabs + Mass*Mass));

    AntiMom->setPx(-px); AntiMom->setPy(-py); AntiMom->setPz(-pz);
    AntiMom->setE (sqrt(Pabs*Pabs + AntiMass*AntiMass));
    }
    
//***********************************************************************************************************    
   
void G4VLongitudinalStringDecay::SetSigmaTransverseMomentum(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetSigmaTransverseMomentum after FragmentString() not allowed");
	} else {
		SigmaQT = aValue;
	}
}

void G4VLongitudinalStringDecay::SetStrangenessSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetStrangenessSuppression after FragmentString() not allowed");
	} else {
		StrangeSuppress = aValue;
	}
}

void G4VLongitudinalStringDecay::SetDiquarkSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetDiquarkSuppression after FragmentString() not allowed");
	} else {
		DiquarkSuppress = aValue;
	}
}

void G4VLongitudinalStringDecay::SetDiquarkBreakProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetDiquarkBreakProbability after FragmentString() not allowed");
	} else {
		DiquarkBreakProb = aValue;
	}
}
		

void G4VLongitudinalStringDecay::SetVectorMesonProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetVectorMesonProbability after FragmentString() not allowed");
	} else {
		pspin_meson = aValue;
	}
}

void G4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetSpinThreeHalfBarionProbability after FragmentString() not allowed");
	} else {
		pspin_barion = aValue;
	}
}

void G4VLongitudinalStringDecay::SetScalarMesonMixings(G4std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetScalarMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      G4Exception("4VLongitudinalStringDecay::SetScalarMesonMixings( argument Vector too small");
	  pmix_meson1[0] = aVector[0];
	  pmix_meson1[1] = aVector[1];
	  pmix_meson1[2] = aVector[2];
	  pmix_meson1[3] = aVector[3];
	  pmix_meson1[4] = aVector[4];
	  pmix_meson1[5] = aVector[5];
	}
}

void G4VLongitudinalStringDecay::SetVectorMesonMixings(G4std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		G4Exception("4VLongitudinalStringDecay::SetVectorMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      G4Exception("4VLongitudinalStringDecay::SetVectorMesonMixings( argument Vector too small");
	  pmix_meson0[0] = aVector[0];
	  pmix_meson0[1] = aVector[1];
	  pmix_meson0[2] = aVector[2];
	  pmix_meson0[3] = aVector[3];
	  pmix_meson0[4] = aVector[4];
	  pmix_meson0[5] = aVector[5];
	}
}	
