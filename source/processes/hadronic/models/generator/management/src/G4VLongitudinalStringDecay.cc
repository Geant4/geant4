// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VLongitudinalStringDecay.cc,v 1.6 1998/12/01 15:41:02 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
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
   
   SigmaQT = 0.5 * GeV;
   
   StrangeSuppress  = 0.44;
   DiquarkSuppress  = 0.1;
   DiquarkBreakProb = 0.1;
   
   SmoothParam      = 0.9; 
   StringLoopInterrupt    = 200;
   ClusterLoopInterrupt   = 500;
   }

G4VLongitudinalStringDecay::~G4VLongitudinalStringDecay()
   {
   }

//********************************************************************************
// Operators

//const G4VLongitudinalStringDecay & G4VLongitudinalStringDecay::operator=(const G4VLongitudinalStringDecay &right)
//    {
//    }

int G4VLongitudinalStringDecay::operator==(const G4VLongitudinalStringDecay &right) const
    {
    return !memcmp(this, &right, sizeof(G4VLongitudinalStringDecay));
    }

int G4VLongitudinalStringDecay::operator!=(const G4VLongitudinalStringDecay &right) const
    {
    return memcmp(this, &right, sizeof(G4VLongitudinalStringDecay));
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
   //... pmix_meson0[] is quark mixing parameters for mesons with spin = 1
    const G4double pmix_meson1[] = {0.5, 0.,   0.5,   0., 1.0, 1.0}; 
   //... pmix_meson1[] is quark mixing parameters for mesons with spin = 3 
    const G4double pmix_meson0[] = {0.5, 0.25, 0.5, 0.25, 1.0, 0.5}; 
   //... pspin_meson is probability to create vector meson 
    const G4double pspin_meson = 0.5;
   //... pspin_barion is probability to create 3/2 barion 
    const G4double pspin_barion = 0.5;

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

      G4int kfld = max(kfla,kflb);
            kfld = max(kfld,kflc);
      G4int kflf = min(kfla,kflb);
            kflf = min(kflf,kflc);

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
   
class SideOfString
   {      
public:   
   G4int    Encoding;
   G4double w;
   G4double Px;
   G4double Py;

   SideOfString(G4int Encoding, G4double w);
   SideOfString(void) {};
   void Init(G4int Encoding, G4double w);
   };

SideOfString::SideOfString(G4int Encoding, G4double w) 
   { 
   this->Encoding = Encoding; 
   this->w = w; 
   Px = Py = 0; 
   }

void SideOfString::Init(G4int Encoding, G4double w) 
   { 
   this->Encoding = Encoding; 
   this->w = w; 
   Px = Py = 0; 
   }

//----------------------------------------------------------------------------------------------------------


G4KineticTrackVector* G4VLongitudinalStringDecay::FragmentString(const G4ExcitedString& theString)
   {
   G4double InitialStringMass = theString.Get4Momentum().mag();
   G4int LeftEncoding  = theString.GetLeftParton()->GetPDGcode();
   G4int RightEncoding = theString.GetRightParton()->GetPDGcode();

   G4double SumMassHeavy;
   G4double SumMassLight;
   G4double ResidualStringMass2 = InitialStringMass*InitialStringMass;
   G4double ResidualStringPx;
   G4double ResidualStringPy;
   G4KineticTrackVector* RightVector = new G4KineticTrackVector; // Create empty output vector
   G4KineticTrackVector* LeftVector  = new G4KineticTrackVector; // Create empty output vector
   G4ParticleDefinition* pHadron;
   G4int cStringLoopInterrupt = 0;

   // Check string decay threshold
       SumMassLight = MassCut;
       if (!(abs(LeftEncoding) > 1000 && abs(RightEncoding) > 1000)) 
          {
          //... string is q --qbar or q--qq type: Build a stable hadron
          // spin 0 meson or spin 1/2 barion will be built

          G4int theSpinLight = (abs(LeftEncoding) < 1000 && abs(RightEncoding) < 1000)? 1: 2;
          pHadron = CreateHadron(LeftEncoding, RightEncoding, true, theSpinLight);
          SumMassLight += pHadron->GetPDGMass();
          if (ResidualStringMass2 <= SumMassLight*SumMassLight)
	      {
	      // Substitute string by light hadron
	      G4ThreeVector Mom3 = theString.Get4Momentum().vect();
	      G4LorentzVector Mom(Mom3, sqrt(Mom3.mag2() + sqr(pHadron->GetPDGMass())));
              delete RightVector;
              LeftVector->insert(new G4KineticTrack(pHadron, 0, theString.GetPosition(), Mom));
	      return LeftVector;
	      }
          }
       else 
	  {
	  //... string is qq--qqbar type: Build two stable hadrons,
	  //... but we need extra uubar or ddbar quark pair
          G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
          if (LeftEncoding < 0) iflc = -iflc;

          //... theSpin = 2; spin 1/2 baryons will be built
          pHadron = CreateHadron(LeftEncoding,  iflc, true, 2);
          SumMassLight += pHadron->GetPDGMass();

          //... theSpin = 2; spin 1/2 baryons will be built
          G4ParticleDefinition* pAntiHadron = CreateHadron(RightEncoding, -iflc, true, 2);
          SumMassLight += pAntiHadron->GetPDGMass();
          if (ResidualStringMass2 <= SumMassLight*SumMassLight)
	      {
	      // Substitute string two light hadron	      
           
	      G4LorentzVector  Mom, AntiMom;
	      Sample4Momentum(&Mom, pHadron->GetPDGMass(), &AntiMom, pAntiHadron->GetPDGMass(), theString.Get4Momentum().mag());
	      LeftVector->insert(new G4KineticTrack(pHadron, 0, theString.GetPosition(), Mom));
	      LeftVector->insert(new G4KineticTrack(pAntiHadron, 0, theString.GetPosition(), AntiMom));
              G4ThreeVector Velocity = theString.Get4Momentum().boostVector();
              delete RightVector;
	      LeftVector->Boost(Velocity);          
	      return LeftVector;
	      }
          }

   SideOfString Left, Right;
   SideOfString* Decay;
   SideOfString* Stable;

LSTART:
   LeftVector->clearAndDestroy();
   RightVector->clearAndDestroy();
   ResidualStringMass2  = InitialStringMass*InitialStringMass;
   Left.Init(LeftEncoding, InitialStringMass);
   Right.Init(RightEncoding, InitialStringMass);
   if (cStringLoopInterrupt++ >= StringLoopInterrupt)
       goto LCLUSTER;
   //.. here starts everything (if the string break up didn't work)
   while(1)
      {     
      SumMassHeavy = MassCut;
      if (!(abs(Left.Encoding) > 1000 && abs(Right.Encoding) > 1000)) 
         {
         //... string is q --qbar or q--qq type: Build a stable hadron
         // spin 1 meson or spin 3/2 barion will be built
         G4int theSpinHeavy = (abs(Left.Encoding) < 1000 && abs(Right.Encoding) < 1000)? 3: 4;
         pHadron = CreateHadron(Left.Encoding, Right.Encoding, true, theSpinHeavy);
         SumMassHeavy += pHadron->GetPDGMass();
         }
      else 
	 {
	 //... string is qq--qqbar type: Build two stable hadrons,
	 //... but we need extra uubar or ddbar quark pair
         G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
         if (Left.Encoding < 0) iflc = -iflc;

         //... theSpin = 4; spin 3/2 baryons will be built
         pHadron = CreateHadron(Left.Encoding,  iflc, true, 4);
         SumMassHeavy += pHadron->GetPDGMass();

         //... theSpin = 4; spin 3/2 baryons will be built
         pHadron = CreateHadron(Right.Encoding, -iflc, true, 4);
         SumMassHeavy += pHadron->GetPDGMass();
         }
      if (ResidualStringMass2 <= SumMassHeavy*SumMassHeavy)
         break; // Goto LastCluster

      //...the main iteration loop for the fragmentation is started here:   
      //... stable and non-stable hadrons can be produced ... update theGivenSpin flag
      //...theGivenSpin = false;
      //...Initialize di_quark_break flag

      Decay  = &Left;
      Stable = &Right;

       //... random choice of string end quark 
       G4int Side = (G4UniformRand() < 0.5)? 1: -1;
       if (Side < 0)
          {
          Decay  = &Right; 
          Stable = &Left;
          }
       //... if string end is  a quark  
       G4int signDecayEncoding = Decay->Encoding/abs(Decay->Encoding);
       G4int QuarkEncoding;
       G4int NewDecayEncoding;
       if (abs(Decay->Encoding) < 1000) 
          {    
          //... it is quark ...  Generate q,qbar or qq,qqbar pair    
          if (G4UniformRand() >= DiquarkSuppress)
             {
             //... q,qbar pair is choosen, sample quark flavor ifln 
             QuarkEncoding    = -signDecayEncoding*SampleQuarkFlavor();
             NewDecayEncoding = -QuarkEncoding;
             }
          else
             {
             //... sample quarks and Construct Diquark ifln 
             G4int i1  = SampleQuarkFlavor();
             G4int i2  = SampleQuarkFlavor();
             G4int i10 = max(i1,i2);
             G4int i20 = min(i1,i2);
             G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
             QuarkEncoding = signDecayEncoding * (i10 * 1000 + i20 * 100 + spin);
             NewDecayEncoding = -QuarkEncoding;
             }    
          } 
       else   
          { 
	  //... it is Diquark ...  
	  //... can Diquark break or not? 
          if (G4UniformRand() < DiquarkBreakProb)
             {
             //... Diquark break

	     G4int StableQuarkEncoding = Decay->Encoding/1000;
	     Decay->Encoding = (Decay->Encoding/100)%10;
	     if (G4UniformRand() < 0.5)
                {
                G4int Swap = StableQuarkEncoding;
                StableQuarkEncoding = Decay->Encoding;
                Decay->Encoding = Swap;
                }
             QuarkEncoding = -Decay->Encoding/abs(Decay->Encoding)*SampleQuarkFlavor();
             //... Build new Diquark
             G4int i10  = max(abs(QuarkEncoding), abs(StableQuarkEncoding));
             G4int i20  = min(abs(QuarkEncoding), abs(StableQuarkEncoding));
             G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
             NewDecayEncoding = -QuarkEncoding/abs(QuarkEncoding)*(i10 * 1000 + i20 * 100 + spin);
             }
          else
             {
             //... Diquark does not break 
             QuarkEncoding = signDecayEncoding*SampleQuarkFlavor();
             NewDecayEncoding = -QuarkEncoding;
             }
          }
       //...  Construct produced hadron: encoding, mass and transverse momentum 

       pHadron = CreateHadron(Decay->Encoding, QuarkEncoding, false,  1);
       Decay->Encoding = NewDecayEncoding;

       G4double HadronMass = pHadron->GetPDGMass();

       // calculate and assign hadron transverse momentum component HadronPx andHadronPy
       G4double thePx, thePy;
       SampleQuarkPt(&thePx, &thePy);
       G4double HadronPx = Decay->Px + thePx;
       G4double HadronPy = Decay->Py + thePy;
       Decay->Px = -thePx;   
       Decay->Py = -thePy;   
       //...  sample z to define hadron longitudinal momentum and energy
       //... but first check the available phase space
       G4double DecayQuarkMass  = FindParticle(Decay->Encoding)->GetPDGMass();
       G4double DecayQuarkMass2 = DecayQuarkMass*DecayQuarkMass;
       G4double HadronMass2T = HadronMass*HadronMass + HadronPx*HadronPx + HadronPy*HadronPy;
       if (DecayQuarkMass2 + HadronMass2T >= SmoothParam*Decay->w*Stable->w) 
          goto LSTART;//continue;  while (2)

       //... then compute allowed z region  z_min <= z <= z_max 
 
       G4double zMin = HadronMass2T/(Decay->w*Stable->w);
       G4double zMax = 1. - DecayQuarkMass2/(Decay->w*Stable->w);
       if (zMin >= zMax) 
	goto LSTART;
       
       G4double z = GetLightConeZ(zMin, zMax, QuarkEncoding, pHadron, HadronPx, HadronPy);      
       
       //... now compute hadron longitudinal momentum and energy
       //longitudinal hadron momentum component HadronPz

       G4double HadronPz = (z * Decay->w - HadronMass2T/(z * Decay->w))*0.5;
       HadronPz *= Side;
       //total hadron energy HadronE
       G4double HadronE  = (z * Decay->w + HadronMass2T/(z * Decay->w))*0.5;

       //...update  (after hadron separation) string light-cone variables
       Left.w  -= HadronE + HadronPz;
       Right.w -= HadronE - HadronPz;

       G4LorentzVector a4Momentum;
       G4ThreeVector   Pos;
       a4Momentum.setPx(HadronPx);a4Momentum.setPy(HadronPy);a4Momentum.setPz(HadronPz);
       a4Momentum.setE (HadronE );
       G4KineticTrack* Hadron = new G4KineticTrack(pHadron, 0, Pos, a4Momentum);
       if (Side > 0)
	  LeftVector->insert(Hadron);
       else
	  RightVector->insert(Hadron); 
       ResidualStringPx  = Left.Px + Right.Px;
       ResidualStringPy  = Left.Py + Right.Py;
       ResidualStringMass2 = Left.w*Right.w - ResidualStringPx*ResidualStringPx- ResidualStringPy*ResidualStringPy;
       SumMassLight = MassCut;
       if (!(abs(Left.Encoding) > 1000 && abs(Right.Encoding) > 1000)) 
          {
          //... string is q --qbar or q--qq type: Build a stable hadron
          // spin 0 meson or spin 1/2 barion will be built

          G4int theSpinLight = (abs(Left.Encoding) < 1000 && abs(Right.Encoding) < 1000)? 1: 2;
          pHadron = CreateHadron(Left.Encoding, Right.Encoding, true, theSpinLight);
          SumMassLight += pHadron->GetPDGMass();
          }
       else 
	  {
	  //... string is qq--qqbar type: Build two stable hadrons,
	  //... but we need extra uubar or ddbar quark pair
          G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
          if (Left.Encoding < 0) iflc = -iflc;

          //... theSpin = 2; spin 1/2 baryons will be built
          pHadron = CreateHadron(Left.Encoding,  iflc, true, 2);
          SumMassLight += pHadron->GetPDGMass();

          //... theSpin = 2; spin 1/2 baryons will be built
          pHadron = CreateHadron(Right.Encoding, -iflc, true, 2);
          SumMassLight += pHadron->GetPDGMass();
          }
       if (ResidualStringMass2 <= SumMassLight*SumMassLight)
	  goto LSTART;
       } 
LCLUSTER:
    //... perform last cluster decay
    ResidualStringPx  = Left.Px + Right.Px;
    ResidualStringPy  = Left.Py + Right.Py;
    G4double ResidualMass     = sqrt(ResidualStringMass2);
    G4double ResidualStringPz = (Left.w - Right.w)/2.;
    G4double ResidualStringE  = (Left.w + Right.w)/2.;
    G4ThreeVector ClusterVel(ResidualStringPx/ResidualStringE, ResidualStringPy/ResidualStringE, ResidualStringPz/ResidualStringE);
    G4ParticleDefinition* pLastHadron;
    Decay  = &Left;
    Stable = &Right;
    G4double ClusterMassCut = ClusterMass;
    G4int cClusterInterrupt = 0;    
    do
       {
       if (cClusterInterrupt++ >= ClusterLoopInterrupt)
          {
          if (cStringLoopInterrupt < StringLoopInterrupt)
             goto LSTART;
          delete RightVector;
          return LeftVector;
          }
       //... if any ifl is a Diquark 
       G4int signDecayEncoding = Decay->Encoding/abs(Decay->Encoding);
       G4int QuarkEncoding;
       if (!(abs(Decay->Encoding) > 1000  || abs(Stable->Encoding) > 1000)) 
          {
	  //... there are quarks on cluster ends
	  //... randomly choose q,qbar pair or qq,qqbar pair
          if (G4UniformRand() > DiquarkSuppress) 
             {
             //... q,qbar pair is choosen, sample quark flavor ifln 
             QuarkEncoding = -signDecayEncoding*SampleQuarkFlavor();
             }    
          else 
             {   
	     //... sample quarks and Construct Diquark ifln 
	     G4int i1 = SampleQuarkFlavor();
	     G4int i2 = SampleQuarkFlavor();
	     G4int i10 = max(i1,i2);
	     G4int i20 = min(i1,i2);
	     G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1: 3;
	     QuarkEncoding = signDecayEncoding*(i10 * 1000 + i20 * 100 + spin);
             }
          }
       else               
          {
	  //... there is a Diquark on cluster ends
	  //...  randomly choose q,qbar pair  
	  QuarkEncoding = ((abs(Decay->Encoding)  < 1000)? -signDecayEncoding : signDecayEncoding)*SampleQuarkFlavor();
          }
       //...  encodings and masses of hadrons  
       //...  theGivenSpin = false; //means any hadron (with given quark content) can be built
       //...  theSpin = 1;
       pHadron = CreateHadron(Decay->Encoding,  QuarkEncoding, false,  1);
       pLastHadron = CreateHadron(Stable->Encoding, -QuarkEncoding, false,  1);

       //... repeat procedure, if mass of cluster is too low to produce hadrons
       //... ClusterMassCut = 0.15*GeV model parameter
       if (QuarkEncoding < 3)
               ClusterMassCut = 0.;
       } 
    while (ResidualMass <= pHadron->GetPDGMass() + pLastHadron->GetPDGMass()  + ClusterMassCut);

    //... compute hadron momenta and energies   
    G4LorentzVector  Mom, AntiMom;
    G4ThreeVector    Pos;
    Sample4Momentum(&Mom, pHadron->GetPDGMass(), &AntiMom, pLastHadron->GetPDGMass(), ResidualMass);
    Mom.boost(ClusterVel);
    AntiMom.boost(ClusterVel);
    LeftVector->insert(new G4KineticTrack(pHadron, 0, Pos, Mom));
    LeftVector->insert(new G4KineticTrack(pLastHadron, 0, Pos, AntiMom));

    // After the next procedure LeftVector will be contained List of Hadrons
    //        in correct order.
    while(!RightVector->isEmpty())
        LeftVector->insert(RightVector->removeLast());
    delete RightVector;

    CalculateHadronTimePosition(InitialStringMass, LeftVector);

    // Make backward lorenz boost 
    G4LorentzRotation toObserverFrame(theString.Get4Momentum().boostVector());
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
    
//*****************************************************************************************************

G4ParticleDefinition* G4VLongitudinalStringDecay::FindParticle(G4int Encoding) 
   {
   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);
   if (ptr == NULL)
       {
       G4cout << "Particle with encoding "<<Encoding<<" does not exist!!!"<<endl;
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
