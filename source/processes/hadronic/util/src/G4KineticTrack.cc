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
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, A. Feliciello, 20th May 1998
// -----------------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
//#include <cmath>

#include "Randomize.hh"
#include "G4SimpleIntegration.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4GeneralPhaseSpaceDecay.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzRotation.hh"
#include "G4SampleResonance.hh"
#include "G4Integrator.hh"
#include "G4KaonZero.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4AntiKaonZero.hh"

#include "G4HadTmpUtil.hh"

//
// Some static clobal for integration
//

static G4ThreadLocal G4double  G4KineticTrack_Gmass, G4KineticTrack_xmass1;

//
//   Default constructor
//

G4KineticTrack::G4KineticTrack() :
                theDefinition(0),
                theFormationTime(0),
                thePosition(0),
                the4Momentum(0),
		theFermi3Momentum(0),
		theTotal4Momentum(0),
		theNucleon(0),
                nChannels(0),
                theActualMass(0),            
                theActualWidth(0),            
                theDaughterMass(0),
                theDaughterWidth(0),
		theStateToNucleus(undefined),
		theProjectilePotential(0),
                theCreatorModel(-1),
                theParentResonanceDef(nullptr),
                theParentResonanceID(0)
{
////////////////
//    DEBUG   //
////////////////

/*
 G4cerr << G4endl << G4endl << G4endl;
 G4cerr << "   G4KineticTrack default constructor invoked! \n";
 G4cerr << "   =========================================== \n" << G4endl;
*/
}



//
//   Copy constructor
//

G4KineticTrack::G4KineticTrack(const G4KineticTrack &right) : G4VKineticNucleon()
{
 theDefinition = right.GetDefinition();
 theFormationTime = right.GetFormationTime();
 thePosition = right.GetPosition();
 the4Momentum = right.GetTrackingMomentum();
 theFermi3Momentum = right.theFermi3Momentum;
 theTotal4Momentum = right.theTotal4Momentum;
 theNucleon=right.theNucleon;
 nChannels = right.GetnChannels();
 theActualMass = right.GetActualMass();
 theActualWidth = new G4double[nChannels];
 for (G4int i = 0; i < nChannels; i++)
  {
    theActualWidth[i] = right.theActualWidth[i];
  }
 theDaughterMass = 0;
 theDaughterWidth = 0;
 theStateToNucleus = right.theStateToNucleus;
 theProjectilePotential = right.theProjectilePotential;
 theCreatorModel = right.GetCreatorModelID();
 theParentResonanceDef = right.GetParentResonanceDef();
 theParentResonanceID = right.GetParentResonanceID();

////////////////
//    DEBUG   //
////////////////

/*
 G4cerr << G4endl << G4endl << G4endl;
 G4cerr << "   G4KineticTrack copy constructor invoked! \n";
 G4cerr << "   ======================================== \n" <<G4endl;
*/
}


//
//   By argument constructor
//

G4KineticTrack::G4KineticTrack(const G4ParticleDefinition* aDefinition,
                               G4double aFormationTime,
                               const G4ThreeVector& aPosition,
                               const G4LorentzVector& a4Momentum) :
                theDefinition(aDefinition),
		theFormationTime(aFormationTime),
                thePosition(aPosition),
                the4Momentum(a4Momentum),
		theFermi3Momentum(0),
		theTotal4Momentum(a4Momentum),
		theNucleon(0),
		theStateToNucleus(undefined),
		theProjectilePotential(0),
                theCreatorModel(-1),
                theParentResonanceDef(nullptr),
                theParentResonanceID(0)
{
  if(G4KaonZero::KaonZero() == theDefinition ||
    G4AntiKaonZero::AntiKaonZero() == theDefinition)
  {
    if(G4UniformRand()<0.5)
    {
      theDefinition = G4KaonZeroShort::KaonZeroShort();
    }
    else
    {
      theDefinition = G4KaonZeroLong::KaonZeroLong();
    }
  }

//
//      Get the number of decay channels
//

 G4DecayTable* theDecayTable = theDefinition->GetDecayTable();
 if (theDecayTable != 0)
    {
     nChannels = theDecayTable->entries();

    }
 else
    {
     nChannels = 0;
    }  

//
//   Get the actual mass value
//

 theActualMass = GetActualMass();

//
//   Create an array to Store the actual partial widths 
//   of the decay channels
//

  theDaughterMass = 0;
  theDaughterWidth = 0;
  theActualWidth = 0;
  G4bool * theDaughterIsShortLived = 0;
  
  if(nChannels!=0) theActualWidth = new G4double[nChannels];

  //  cout << " ****CONSTR*** ActualMass ******* " << theActualMass << G4endl;
  G4int index;
  for (index = nChannels - 1; index >= 0; --index)
    {
      G4VDecayChannel* theChannel = theDecayTable->GetDecayChannel(index);
      G4int nDaughters = theChannel->GetNumberOfDaughters();
      G4double theMotherWidth;
      if (nDaughters == 2 || nDaughters == 3) 
	{
          G4double thePoleMass  = theDefinition->GetPDGMass();
          theMotherWidth = theDefinition->GetPDGWidth();
          G4double thePoleWidth = theChannel->GetBR()*theMotherWidth;
          const G4ParticleDefinition* aDaughter;
          theDaughterMass = new G4double[nDaughters];
          theDaughterWidth = new G4double[nDaughters];
	  theDaughterIsShortLived = new G4bool[nDaughters];
          for (G4int n = 0; n < nDaughters; ++n)
	    {
              aDaughter = theChannel->GetDaughter(n);
              theDaughterMass[n] = aDaughter->GetPDGMass();
              theDaughterWidth[n] = aDaughter->GetPDGWidth();
	      theDaughterIsShortLived[n] = aDaughter->IsShortLived();
	    }     
	  
//
//           Check whether both the decay products are stable
//

          G4double theActualMom = 0.0;
          G4double thePoleMom = 0.0;
	  G4SampleResonance aSampler;
	  if (nDaughters==2) 
	    {
	      if ( !theDaughterIsShortLived[0] && !theDaughterIsShortLived[1] )
		{
		  
		  //              G4cout << G4endl << "Both the " << nDaughters <<
		  //                              " decay products are stable!";
		  //	       	       cout << " LB: Both decay products STABLE !" << G4endl;
		  //	       	       cout << " parent:     " << theChannel->GetParentName() << G4endl;
		  //	       	       cout << " particle1:  " << theChannel->GetDaughterName(0) << G4endl;
		  //	       	       cout << " particle2:  " << theChannel->GetDaughterName(1) << G4endl;
		  
		  theActualMom = EvaluateCMMomentum(theActualMass, 
						    theDaughterMass);   
		  thePoleMom = EvaluateCMMomentum(thePoleMass, 
						  theDaughterMass);
		  //	      cout << G4endl;
		  //	      cout << " LB: ActualMass/DaughterMass  " << theActualMass << "   " << theDaughterMass << G4endl; 
		  //	      cout << " LB: ActualMom " << theActualMom << G4endl;
		  //	      cout << " LB: PoleMom   " << thePoleMom << G4endl;
		  //	      cout << G4endl;
		}
	      else if ( !theDaughterIsShortLived[0] && theDaughterIsShortLived[1] )   
		{
		  
		  //              G4cout << G4endl << "Only the first of the " << nDaughters <<" decay products is stable!";
		  //	       cout << " LB: only the first decay product is STABLE !" << G4endl;
		  //	       cout << " parent:     " << theChannel->GetParentName() << G4endl;
		  //	       cout << " particle1:  " << theChannel->GetDaughterName(0) << G4endl;
		  //	       cout << " particle2:  " << theChannel->GetDaughterName(1) << G4endl;
		  
// global variable definition
		  G4double lowerLimit = aSampler.GetMinimumMass(theChannel->GetDaughter(1));
		  theActualMom = IntegrateCMMomentum(lowerLimit);
		  thePoleMom = IntegrateCMMomentum(lowerLimit, thePoleMass);
		  //	      cout << " LB Parent Mass = " <<  G4KineticTrack_Gmass << G4endl;
		  //	      cout << " LB Actual Mass = " << theActualMass << G4endl;
		  //	      cout << " LB Daughter1 Mass = " <<  G4KineticTrack_Gmass1 << G4endl;
		  //	      cout << " LB Daughter2 Mass = " <<  G4KineticTrack_Gmass2 << G4endl;
		  //	      cout << " The Actual Momentum = " << theActualMom << G4endl;
		  //	      cout << " The Pole Momentum   = " << thePoleMom << G4endl;
		  //	      cout << G4endl;
		  
		}        
	      else if ( theDaughterIsShortLived[0] && !theDaughterIsShortLived[1] )   
		{
		  
		  //              G4cout << G4endl << "Only the second of the " << nDaughters <<
		  //                              " decay products is stable!";
		  //	       	       cout << " LB: only the second decay product is STABLE !" << G4endl;
		  //	       cout << " parent:     " << theChannel->GetParentName() << G4endl;
		  //	       cout << " particle1:  " << theChannel->GetDaughterName(0) << G4endl;
		  //	       cout << " particle2:  " << theChannel->GetDaughterName(1) << G4endl;
		  
//
//               Swap the content of the theDaughterMass and theDaughterWidth arrays!!!
//
		  
		  G4SwapObj(theDaughterMass, theDaughterMass + 1);
		  G4SwapObj(theDaughterWidth, theDaughterWidth + 1);
		  
// global variable definition
		  G4double lowerLimit = aSampler.GetMinimumMass(theChannel->GetDaughter(0));
		  theActualMom = IntegrateCMMomentum(lowerLimit);
		  thePoleMom = IntegrateCMMomentum(lowerLimit, thePoleMass);
		  //	      cout << " LB Parent Mass = " <<  G4KineticTrack_Gmass << G4endl;
		  //	      cout << " LB Actual Mass = " << theActualMass << G4endl;
		  //	      cout << " LB Daughter1 Mass = " <<  G4KineticTrack_Gmass1 << G4endl;
		  //	      cout << " LB Daughter2 Mass = " <<  G4KineticTrack_Gmass2 << G4endl;
		  //	      cout << " The Actual Momentum = " << theActualMom << G4endl;
		  //	      cout << " The Pole Momentum   = " << thePoleMom << G4endl;
		  //              cout << G4endl;
		  
		}        
	      else if ( theDaughterIsShortLived[0] && theDaughterIsShortLived[1] )   
		{
	       
//              G4cout << G4endl << "Both the " << nDaughters <<
//                              " decay products are resonances!";
		  //	       cout << " LB: both decay products are RESONANCES !" << G4endl;
		  //	       cout << " parent:     " << theChannel->GetParentName() << G4endl;
		  //	       cout << " particle1:  " << theChannel->GetDaughterName(0) << G4endl;
		  //	       cout << " particle2:  " << theChannel->GetDaughterName(1) << G4endl;
		  
// global variable definition
		  G4KineticTrack_Gmass = theActualMass;
		  theActualMom = IntegrateCMMomentum2();
		  G4KineticTrack_Gmass = thePoleMass;
		  thePoleMom = IntegrateCMMomentum2();
		  //	      cout << " LB Parent Mass = " <<  G4KineticTrack_Gmass << G4endl;
		  //	      cout << " LB Daughter1 Mass = " <<  G4KineticTrack_Gmass1 << G4endl;
		  //	      cout << " LB Daughter2 Mass = " <<  G4KineticTrack_Gmass2 << G4endl;
		  //              cout << " The Actual Momentum = " << theActualMom << G4endl;
		  //              cout << " The Pole Momentum   = " << thePoleMom << G4endl;
		  //              cout << G4endl;
		  
		}        
	    } 
	  else  // (nDaughter==3)
	    {
	      
	      G4int nShortLived = 0;
	      if ( theDaughterIsShortLived[0] ) 
		{ 
		  ++nShortLived; 
		}
	      if ( theDaughterIsShortLived[1] )
		{ 
		  ++nShortLived; 
		  G4SwapObj(theDaughterMass, theDaughterMass + 1);
		  G4SwapObj(theDaughterWidth, theDaughterWidth + 1);	    
		}
	      if ( theDaughterIsShortLived[2] )
		{ 
		  ++nShortLived; 
		  G4SwapObj(theDaughterMass, theDaughterMass + 2);
		  G4SwapObj(theDaughterWidth, theDaughterWidth + 2);	    
		}
	      if ( nShortLived == 0 ) 
		{
		  theDaughterMass[1]+=theDaughterMass[2];
		  theActualMom = EvaluateCMMomentum(theActualMass, 
						    theDaughterMass);   
		  thePoleMom = EvaluateCMMomentum(thePoleMass, 
						  theDaughterMass);
		}
//	      else if ( nShortLived == 1 )
	      else if ( nShortLived >= 1 )
		{ 
		  // need the shortlived particle in slot 1! (very bad style...)
		  G4SwapObj(theDaughterMass, theDaughterMass + 1);
		  G4SwapObj(theDaughterWidth, theDaughterWidth + 1);	    
		  theDaughterMass[0] += theDaughterMass[2];
		  theActualMom = IntegrateCMMomentum(0.0);
		  thePoleMom = IntegrateCMMomentum(0.0, thePoleMass);
		}
//	      else
//		{
//		  throw G4HadronicException(__FILE__, __LINE__,  ("can't handle more than one shortlived in 3 particle output channel");
//		}     
	      
	    }
	  
	  //if(nDaughters<3) theChannel->GetAngularMomentum(); 
	  G4double theMassRatio = thePoleMass / theActualMass;
          G4double theMomRatio = theActualMom / thePoleMom;
	  // VI 11.06.2015: for l=0 one not need use pow
          //G4double l=0;
          //theActualWidth[index] = thePoleWidth * theMassRatio *
          //                        std::pow(theMomRatio, (2 * l + 1)) *
          //                        (1.2 / (1+ 0.2*std::pow(theMomRatio, (2 * l))));
          theActualWidth[index] = thePoleWidth * theMassRatio *
                                  theMomRatio;
          delete [] theDaughterMass;
	  theDaughterMass = 0;
          delete [] theDaughterWidth;
	  theDaughterWidth = 0;
	  delete [] theDaughterIsShortLived;
          theDaughterIsShortLived = 0;
	}
      
      else //  nDaughter = 1 ( e.g. K0  decays 50% to Kshort, 50% Klong 
	{
	  theMotherWidth = theDefinition->GetPDGWidth();
	  theActualWidth[index] = theChannel->GetBR()*theMotherWidth;
	}
    }

////////////////
//    DEBUG   //
////////////////

// for (G4int y = nChannels - 1; y >= 0; --y)
//     {
//      G4cout << G4endl << theActualWidth[y];
//     }
// G4cout << G4endl << G4endl << G4endl;

 /*
 G4cerr << G4endl << G4endl << G4endl;
 G4cerr << "   G4KineticTrack by argument constructor invoked! \n";
 G4cerr << "   =============================================== \n" << G4endl;
 */

}

G4KineticTrack::G4KineticTrack(G4Nucleon * nucleon,
			       const G4ThreeVector& aPosition,
                               const G4LorentzVector& a4Momentum)
  :     theDefinition(nucleon->GetDefinition()),
	theFormationTime(0),
	thePosition(aPosition),
	the4Momentum(a4Momentum),
	theFermi3Momentum(nucleon->GetMomentum()),
        theNucleon(nucleon),
	nChannels(0),
	theActualMass(nucleon->GetDefinition()->GetPDGMass()),
	theActualWidth(0),
	theDaughterMass(0),
	theDaughterWidth(0),
	theStateToNucleus(undefined),
	theProjectilePotential(0),
        theCreatorModel(-1),
        theParentResonanceDef(nullptr),
        theParentResonanceID(0) 
{
	theFermi3Momentum.setE(0);
	Set4Momentum(a4Momentum);
}


G4KineticTrack::~G4KineticTrack()
{
 if (theActualWidth != 0) delete [] theActualWidth;
 if (theDaughterMass != 0) delete [] theDaughterMass;
 if (theDaughterWidth != 0) delete [] theDaughterWidth;
}



G4KineticTrack& G4KineticTrack::operator=(const G4KineticTrack& right)
{
 if (this != &right)
    {
     theDefinition = right.GetDefinition();
     theFormationTime = right.GetFormationTime();
     the4Momentum = right.the4Momentum;  
     the4Momentum = right.GetTrackingMomentum();
     theFermi3Momentum = right.theFermi3Momentum;
     theTotal4Momentum = right.theTotal4Momentum;
     theNucleon = right.theNucleon;
     theStateToNucleus = right.theStateToNucleus;
     if (theActualWidth != 0) delete [] theActualWidth;
     nChannels = right.GetnChannels();      
     theActualWidth = new G4double[nChannels];
     for (G4int i = 0; i < nChannels; ++i) theActualWidth[i] = right.theActualWidth[i];
     theCreatorModel = right.GetCreatorModelID();
     theParentResonanceDef = right.GetParentResonanceDef();
     theParentResonanceID = right.GetParentResonanceID();
    }
 return *this;
}



G4bool G4KineticTrack::operator==(const G4KineticTrack& right) const
{
 return (this == & right);
}



G4bool G4KineticTrack::operator!=(const G4KineticTrack& right) const
{
 return (this != & right);
}



G4KineticTrackVector* G4KineticTrack::Decay()
{
//
//   Select a possible decay channel
//
/*
    G4int index1;
    for (index1 = nChannels - 1; index1 >= 0; --index1)
      G4cout << "DECAY Actual Width IND/ActualW " << index1 << "  " << theActualWidth[index1] << G4endl;
      G4cout << "DECAY Actual Mass " << theActualMass << G4endl;
*/
  const G4ParticleDefinition* thisDefinition = this->GetDefinition();
  if(!thisDefinition)
  {
    G4cerr << "Error condition encountered in G4KineticTrack::Decay()"<<G4endl;
    G4cerr << "  track has no particle definition associated."<<G4endl;
    return 0;
  }
  G4DecayTable* theDecayTable = thisDefinition->GetDecayTable();
  if(!theDecayTable)
  {
    G4cerr << "Error condition encountered in G4KineticTrack::Decay()"<<G4endl;
    G4cerr << "  particle definition has no decay table associated."<<G4endl;
    G4cerr << "  particle was "<<thisDefinition->GetParticleName()<<G4endl;
    return 0;
  }
 
 G4int chargeBalance = G4lrint(theDefinition->GetPDGCharge() );     
 G4int baryonBalance = G4lrint(theDefinition->GetBaryonNumber() );
 G4LorentzVector energyMomentumBalance(Get4Momentum());
 G4double theTotalActualWidth = this->EvaluateTotalActualWidth();
 if (theTotalActualWidth !=0)
    {

     //AR-16Aug2016 : Repeat the sampling of the decay channel until is 
     //               kinematically above threshold or a max number of attempts is reached
     G4bool isChannelBelowThreshold = true;
     const G4int maxNumberOfLoops = 10000;
     G4int loopCounter = 0;

     G4int chosench;
     G4String theParentName;
     G4double theParentMass;
     G4double theBR;
     G4int theNumberOfDaughters;
     G4String theDaughtersName1;
     G4String theDaughtersName2;
     G4String theDaughtersName3;
     G4String theDaughtersName4;
     G4double masses[4]={0.,0.,0.,0.};

     do {       

       G4double theSumActualWidth = 0.0;
       G4double* theCumActualWidth = new G4double[nChannels]{};
       for (G4int index = nChannels - 1; index >= 0; --index)
          {
           theSumActualWidth += theActualWidth[index];
           theCumActualWidth[index] = theSumActualWidth;
	   //	 cout << "DECAY Cum. Width " << index << "  " << theCumActualWidth[index] << G4endl;
	  }
       //	 cout << "DECAY Total Width " << theSumActualWidth << G4endl;
       //	 cout << "DECAY Total Width " << theTotalActualWidth << G4endl;
       G4double r = theTotalActualWidth * G4UniformRand();
       G4VDecayChannel* theDecayChannel(0);
       chosench=-1;
       for (G4int index = nChannels - 1; index >= 0; --index)
          {
           if (r < theCumActualWidth[index])
              {
               theDecayChannel = theDecayTable->GetDecayChannel(index);
	       //	     cout << "DECAY SELECTED CHANNEL" << index << G4endl;
               chosench=index;
               break; 
              }
          }

       delete [] theCumActualWidth;
   
       if(!theDecayChannel)
       {
         G4cerr << "Error condition encountered in G4KineticTrack::Decay()"<<G4endl;
         G4cerr << "  decay channel has 0x0 channel associated."<<G4endl;
         G4cerr << "  particle was "<<thisDefinition->GetParticleName()<<G4endl;
         G4cerr << "  channel index "<< chosench << "of "<<nChannels<<"channels"<<G4endl;
         return 0;
       }
       theParentName = theDecayChannel->GetParentName();
       theParentMass = this->GetActualMass();
       theBR = theActualWidth[chosench];
       //     cout << "**BR*** DECAYNEW  " << theBR << G4endl;
       theNumberOfDaughters = theDecayChannel->GetNumberOfDaughters();
       theDaughtersName1 = "";
       theDaughtersName2 = "";
       theDaughtersName3 = "";
       theDaughtersName4 = "";

       for (G4int i=0; i < 4; ++i) masses[i]=0.;
       G4int shortlivedDaughters[4];
       G4int numberOfShortliveds(0);
       G4double SumLongLivedMass(0);
       for (G4int aD=0; aD < theNumberOfDaughters ; ++aD)
       {
          const G4ParticleDefinition* aDaughter = theDecayChannel->GetDaughter(aD);
          masses[aD] = aDaughter->GetPDGMass();
          if ( aDaughter->IsShortLived() ) 
	  {
	      shortlivedDaughters[numberOfShortliveds]=aD;
	      ++numberOfShortliveds;
	  } else {
	      SumLongLivedMass += aDaughter->GetPDGMass();
	  }
		
       }    
       switch (theNumberOfDaughters)
          {
           case 0:
              break;
           case 1:
              theDaughtersName1 = theDecayChannel->GetDaughterName(0);
              theDaughtersName2 = "";
              theDaughtersName3 = "";
              theDaughtersName4 = "";
              break;
           case 2:    
              theDaughtersName1 = theDecayChannel->GetDaughterName(0);	    
              theDaughtersName2 = theDecayChannel->GetDaughterName(1);
              theDaughtersName3 = "";
              theDaughtersName4 = "";
	      if (  numberOfShortliveds == 1) 
	      {   G4SampleResonance aSampler;
                  G4double massmax=theParentMass - SumLongLivedMass;
		  const G4ParticleDefinition * aDaughter=theDecayChannel->GetDaughter(shortlivedDaughters[0]);
	          masses[shortlivedDaughters[0]]= aSampler.SampleMass(aDaughter,massmax);
	      } else if (  numberOfShortliveds == 2) {
	          // choose masses one after the other, start with randomly choosen
	          G4int zero= (G4UniformRand() > 0.5) ? 0 : 1;
		  G4int one = 1-zero;
		  G4SampleResonance aSampler;
		  G4double massmax=theParentMass - aSampler.GetMinimumMass(theDecayChannel->GetDaughter(shortlivedDaughters[one]));
		  const G4ParticleDefinition * aDaughter=theDecayChannel->GetDaughter(shortlivedDaughters[zero]);
		  masses[shortlivedDaughters[zero]]=aSampler.SampleMass(aDaughter,massmax);
		  massmax=theParentMass - masses[shortlivedDaughters[zero]];
		  aDaughter=theDecayChannel->GetDaughter(shortlivedDaughters[one]);
		  masses[shortlivedDaughters[one]]=aSampler.SampleMass(aDaughter,massmax);
	      }
	      break;		
	   case 3:    
              theDaughtersName1 = theDecayChannel->GetDaughterName(0);
              theDaughtersName2 = theDecayChannel->GetDaughterName(1);
              theDaughtersName3 = theDecayChannel->GetDaughterName(2);
              theDaughtersName4 = "";
	      if (  numberOfShortliveds == 1) 
	      {   G4SampleResonance aSampler;
                  G4double massmax=theParentMass - SumLongLivedMass;
	    	  const G4ParticleDefinition * aDaughter=theDecayChannel->GetDaughter(shortlivedDaughters[0]);
	          masses[shortlivedDaughters[0]]= aSampler.SampleMass(aDaughter,massmax);
	      }
	      break;
	   default:    
              theDaughtersName1 = theDecayChannel->GetDaughterName(0);
              theDaughtersName2 = theDecayChannel->GetDaughterName(1);
              theDaughtersName3 = theDecayChannel->GetDaughterName(2);
              theDaughtersName4 = theDecayChannel->GetDaughterName(3);
	      if (  numberOfShortliveds == 1) 
	      {   G4SampleResonance aSampler;
                  G4double massmax=theParentMass - SumLongLivedMass;
	    	  const G4ParticleDefinition * aDaughter=theDecayChannel->GetDaughter(shortlivedDaughters[0]);
	          masses[shortlivedDaughters[0]]= aSampler.SampleMass(aDaughter,massmax);
	      }
              if ( theNumberOfDaughters > 4 ) {
                G4ExceptionDescription ed;
                ed << "More than 4 decay daughters: kept only the first 4" << G4endl;
                G4Exception( "G4KineticTrack::Decay()", "KINTRK5", JustWarning, ed );
              }
	      break;
	  }

          //AR-16Aug2016 : Check whether the sum of the masses of the daughters is smaller than the parent mass.
          //               If this is still not the case, but the max number of attempts has been reached,
          //               then the subsequent call thePhaseSpaceDecayChannel.DecayIt() will throw an exception.
          G4double sumDaughterMasses = 0.0;
          for (G4int i=0; i < 4; ++i) sumDaughterMasses += masses[i];
          if ( theParentMass - sumDaughterMasses > 0.0 ) isChannelBelowThreshold = false;

     } while ( isChannelBelowThreshold && ++loopCounter < maxNumberOfLoops );   /* Loop checking, 16.08.2016, A.Ribon */

//	
//      Get the decay products List
//
     
     G4GeneralPhaseSpaceDecay thePhaseSpaceDecayChannel(theParentName,
                                                        theParentMass,
                                                        theBR,
                                                        theNumberOfDaughters,
                                                        theDaughtersName1,                  
		                                        theDaughtersName2,
		                                        theDaughtersName3,
		                                        theDaughtersName4,
							masses);
     G4DecayProducts* theDecayProducts = thePhaseSpaceDecayChannel.DecayIt();
     if(!theDecayProducts)
     {
       G4ExceptionDescription ed;
       ed << "Error condition encountered: phase-space decay failed." << G4endl
          << "\t the decaying particle is: " << thisDefinition->GetParticleName() << G4endl
          << "\t the channel index is: "<< chosench << " of "<< nChannels << "channels" << G4endl
          << "\t " << theNumberOfDaughters << " daughter particles: "
          << theDaughtersName1 << " " << theDaughtersName2 << " " << theDaughtersName3 << " " 
          << theDaughtersName4 << G4endl;
       G4Exception( "G4KineticTrack::Decay ", "HAD_KINTRACK_001", JustWarning, ed );
       return 0;
     }
		                                        
//
//      Create the kinetic track List associated to the decay products
//      
//      For the decay products of hadronic resonances, we assign as creator model ID
//      the same as their parent
     G4LorentzRotation toMoving(Get4Momentum().boostVector());
     G4DynamicParticle* theDynamicParticle;
     G4double formationTime = 0.0;
     G4ThreeVector position = this->GetPosition();
     G4LorentzVector momentum;
     G4LorentzVector momentumBalanceCMS(0);
     G4KineticTrackVector* theDecayProductList = new G4KineticTrackVector;
     G4int dEntries = theDecayProducts->entries();
     const G4ParticleDefinition * aProduct = 0;
     // Use the integer round mass in keV to get an unique ID for the parent resonance
     G4int uniqueID = static_cast< G4int >( round( Get4Momentum().mag() / CLHEP::keV ) );
     for (G4int i=dEntries; i > 0; --i)
        {
	 theDynamicParticle = theDecayProducts->PopProducts();
         aProduct = theDynamicParticle->GetDefinition();
         chargeBalance -= G4lrint(aProduct->GetPDGCharge() );
         baryonBalance -= G4lrint(aProduct->GetBaryonNumber() );
	 momentumBalanceCMS += theDynamicParticle->Get4Momentum();
         momentum = toMoving*theDynamicParticle->Get4Momentum();
         energyMomentumBalance -= momentum;
         G4KineticTrack* aDaughter = new G4KineticTrack (aProduct,
                                                         formationTime,
                                                         position,
                                                         momentum);
         if (aDaughter != nullptr) 
	   {
             aDaughter->SetCreatorModelID(GetCreatorModelID());
             aDaughter->SetParentResonanceDef(GetDefinition());
             aDaughter->SetParentResonanceID(uniqueID);
           }
         theDecayProductList->push_back(aDaughter);
         delete theDynamicParticle;
        }
     delete theDecayProducts;
     if(std::getenv("DecayEnergyBalanceCheck"))
       std::cout << "DEBUGGING energy balance in cms and lab, charge baryon balance : "
       	         << momentumBalanceCMS << " " 
  	         <<energyMomentumBalance << " " 
  	         <<chargeBalance<<" "
	         <<baryonBalance<<" "
  	         <<G4endl;
     return theDecayProductList;
    }
 else
    {
     return 0;
    }
}

G4double G4KineticTrack::IntegrandFunction1(G4double xmass) const 
{
  G4double mass = theActualMass;   /* the actual mass value */
  G4double mass1 = theDaughterMass[0];
  G4double mass2 = theDaughterMass[1];
  G4double gamma2 = theDaughterWidth[1];
  
  G4double result = (1. / (2 * mass)) *
    std::sqrt(std::max((((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
	     ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))),0.0)) *
    BrWig(gamma2, mass2, xmass);
  return result;
}

G4double G4KineticTrack::IntegrandFunction2(G4double xmass) const
{
  G4double mass = theDefinition->GetPDGMass();   /* the pole mass value */
  G4double mass1 = theDaughterMass[0];
  G4double mass2 = theDaughterMass[1];
  G4double gamma2 = theDaughterWidth[1];
  G4double result = (1. / (2 * mass)) *
    std::sqrt(std::max((((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
	     ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))),0.0)) *
    BrWig(gamma2, mass2, xmass);
 return result;
}

G4double G4KineticTrack::IntegrandFunction3(G4double xmass) const
{
  const G4double mass =  G4KineticTrack_Gmass;   /* the actual mass value */
//  const G4double mass1 = theDaughterMass[0];
  const G4double mass2 = theDaughterMass[1];
  const G4double gamma2 = theDaughterWidth[1];

  const G4double result = (1. / (2 * mass)) *
    std::sqrt(((mass * mass) - (G4KineticTrack_xmass1 + xmass) * (G4KineticTrack_xmass1 + xmass)) *
	 ((mass * mass) - (G4KineticTrack_xmass1 - xmass) * (G4KineticTrack_xmass1 - xmass))) *
    BrWig(gamma2, mass2, xmass);
  return result;
}

G4double G4KineticTrack::IntegrandFunction4(G4double xmass) const
{
  const G4double mass =  G4KineticTrack_Gmass;
  const G4double mass1 = theDaughterMass[0];
  const G4double gamma1 = theDaughterWidth[0];
//  const G4double mass2 = theDaughterMass[1];
  
  G4KineticTrack_xmass1 = xmass;
  
  const G4double theLowerLimit = 0.0;
  const G4double theUpperLimit = mass - xmass;
  const G4int nIterations = 100;
  
  G4Integrator<const G4KineticTrack, G4double(G4KineticTrack::*)(G4double) const> integral;
  G4double result = BrWig(gamma1, mass1, xmass)*
    integral.Simpson(this, &G4KineticTrack::IntegrandFunction3, theLowerLimit, theUpperLimit, nIterations);

  return result;
}

G4double G4KineticTrack::IntegrateCMMomentum(const G4double theLowerLimit) const
{
  const G4double theUpperLimit = theActualMass - theDaughterMass[0];
  const G4int nIterations = 100;
 
  if (theLowerLimit>=theUpperLimit) return 0.0;

  G4Integrator<const G4KineticTrack, G4double(G4KineticTrack::*)(G4double) const> integral;
  G4double theIntegralOverMass2 = integral.Simpson(this, &G4KineticTrack::IntegrandFunction1, 
						   theLowerLimit, theUpperLimit, nIterations);
  return theIntegralOverMass2;
}

G4double G4KineticTrack::IntegrateCMMomentum(const G4double theLowerLimit, const G4double poleMass) const
{
  const G4double theUpperLimit = poleMass - theDaughterMass[0];
  const G4int nIterations = 100;
  
  if (theLowerLimit>=theUpperLimit) return 0.0;

  G4Integrator<const G4KineticTrack, G4double(G4KineticTrack::*)(G4double) const> integral;
  const G4double theIntegralOverMass2 = integral.Simpson(this, &G4KineticTrack::IntegrandFunction2,
						    theLowerLimit, theUpperLimit, nIterations);
  return theIntegralOverMass2;
}


G4double G4KineticTrack::IntegrateCMMomentum2() const
{
  const G4double theLowerLimit = 0.0;
  const G4double theUpperLimit = theActualMass;
  const G4int nIterations = 100;
  
  if (theLowerLimit>=theUpperLimit) return 0.0;
  
  G4Integrator<const G4KineticTrack, G4double(G4KineticTrack::*)(G4double) const> integral;
  G4double theIntegralOverMass2 = integral.Simpson(this, &G4KineticTrack::IntegrandFunction4,
						   theLowerLimit, theUpperLimit, nIterations);
  return theIntegralOverMass2;
}








