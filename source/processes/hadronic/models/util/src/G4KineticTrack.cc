//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include <math.h>

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

//
// Some static clobal for integration
//

static G4double  G4KineticTrack_Gmass, G4KineticTrack_xmass1;

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
		theStateToNucleus(undefined)
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
 G4int i;
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
 for (i = 0; i < nChannels; i++)
  {
    theActualWidth[i] = right.theActualWidth[i];
  }
  theDaughterMass = 0;
  theDaughterWidth = 0;
  theStateToNucleus=right.theStateToNucleus;
 
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

G4KineticTrack::G4KineticTrack(G4ParticleDefinition* aDefinition,
                               G4double aFormationTime,
                               G4ThreeVector aPosition,
                               G4LorentzVector& a4Momentum) :
                theDefinition(aDefinition),
		theFormationTime(aFormationTime),
                thePosition(aPosition),
                the4Momentum(a4Momentum),
		theFermi3Momentum(0),
		theTotal4Momentum(a4Momentum),
		theNucleon(0),
		theStateToNucleus(undefined)
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
  for (index = nChannels - 1; index >= 0; index--)
    {
      G4VDecayChannel* theChannel = theDecayTable->GetDecayChannel(index);
      G4int nDaughters = theChannel->GetNumberOfDaughters();
      G4double theMotherWidth;
      if (nDaughters == 2 || nDaughters == 3) 
	{
          G4double thePoleMass  = theDefinition->GetPDGMass();
          theMotherWidth = theDefinition->GetPDGWidth();
          G4double thePoleWidth = theChannel->GetBR()*theMotherWidth;
          G4ParticleDefinition* aDaughter;
          theDaughterMass = new G4double[nDaughters];
          theDaughterWidth = new G4double[nDaughters];
	  theDaughterIsShortLived = new G4bool[nDaughters];
          G4int n;
          for (n = 0; n < nDaughters; n++)
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
	      
	      int nShortLived = 0;
	      if ( theDaughterIsShortLived[0] ) 
		{ 
		  nShortLived++; 
		}
	      if ( theDaughterIsShortLived[1] )
		{ 
		  nShortLived++; 
		  G4SwapObj(theDaughterMass, theDaughterMass + 1);
		  G4SwapObj(theDaughterWidth, theDaughterWidth + 1);	    
		}
	      if ( theDaughterIsShortLived[2] )
		{ 
		  nShortLived++; 
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
//		  G4Exception ("can't handle more than one shortlived in 3 particle output channel");
//		}     
	      
	    }
	  
          G4double l=0;
	  //if(nDaughters<3) theChannel->GetAngularMomentum(); 
	  G4double theMassRatio = thePoleMass / theActualMass;
          G4double theMomRatio = theActualMom / thePoleMom;
          theActualWidth[index] = thePoleWidth * theMassRatio *
                                  pow(theMomRatio, (2 * l + 1)) *
                                  (1.2 / (1+ 0.2*pow(theMomRatio, (2 * l))));
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

// for (G4int y = nChannels - 1; y >= 0; y--)
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
                                G4ThreeVector aPosition,
                                G4LorentzVector& a4Momentum)
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
	theStateToNucleus(undefined)
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



const G4KineticTrack& G4KineticTrack::operator=(const G4KineticTrack& right)
{
 G4int i;
 if (this != &right)
    {
     theDefinition = right.GetDefinition();
     theFormationTime = right.GetFormationTime();
//     thePosition = right.GetPosition();
     the4Momentum = right.the4Momentum;  
     the4Momentum = right.GetTrackingMomentum();
     theFermi3Momentum = right.theFermi3Momentum;
     theTotal4Momentum = right.theTotal4Momentum;
     theNucleon=right.theNucleon;
     theStateToNucleus=right.theStateToNucleus;
     if (theActualWidth != 0) delete [] theActualWidth;
     nChannels = right.GetnChannels();      
     theActualWidth = new G4double[nChannels];
     for (i = 0; i < nChannels; i++) 
        {
         theActualWidth[i] = right.theActualWidth[i];
        }
    }
 return *this;
}



G4int G4KineticTrack::operator==(const G4KineticTrack& right) const
{
 return (this == & right);
}



G4int G4KineticTrack::operator!=(const G4KineticTrack& right) const
{
 return (this != & right);
}



G4KineticTrackVector* G4KineticTrack::Decay()
{
//
//   Select a possible decay channel
//
  //  G4int index1;
  //  for (index1 = nChannels - 1; index1 >= 0; index1--)
    //  cout << "DECAY Actual Width IND/ActualW " << index1 << "  " << theActualWidth[index1] << G4endl;
    //  cout << "DECAY Actual Mass " << theActualMass << G4endl;
 
 G4double theTotalActualWidth = this->EvaluateTotalActualWidth();
 if (theTotalActualWidth !=0)
    {
     G4int index;
     G4double theSumActualWidth = 0.0;
     G4double* theCumActualWidth = new G4double[nChannels];
     for (index = nChannels - 1; index >= 0; index--)
        {
         theSumActualWidth += theActualWidth[index];
         theCumActualWidth[index] = theSumActualWidth;
	 //	 cout << "DECAY Cum. Width " << index << "  " << theCumActualWidth[index] << G4endl;
	}
     //	 cout << "DECAY Total Width " << theSumActualWidth << G4endl;
     //	 cout << "DECAY Total Width " << theTotalActualWidth << G4endl;
     G4double r = theTotalActualWidth * G4UniformRand();
     G4ParticleDefinition* theDefinition = this->GetDefinition();
     G4DecayTable* theDecayTable = theDefinition->GetDecayTable();
     G4VDecayChannel* theDecayChannel=NULL;
     for (index = nChannels - 1; index >= 0; index--)
        {
         if (r < theCumActualWidth[index])
            {
             theDecayChannel = theDecayTable->GetDecayChannel(index);
	     //	     cout << "DECAY SELECTED CHANNEL" << index << G4endl;
             chosench=index;
             break; 
            }
        }
        
     G4String theParentName = theDecayChannel->GetParentName();
     G4double theParentMass = this->GetActualMass();
     G4double theBR = theActualWidth[index];
     //     cout << "**BR*** DECAYNEW  " << theBR << G4endl;
     //     cout << "**PMass*** DECAYNEW  " << theParentMass << G4endl;
     G4int theNumberOfDaughters = theDecayChannel->GetNumberOfDaughters();
     G4String theDaughtersName1 = "";
     G4String theDaughtersName2 = "";
     G4String theDaughtersName3 = "";     
     switch (theNumberOfDaughters)
        {
         case 0:
            break;
         case 1:
            theDaughtersName1 = theDecayChannel->GetDaughterName(0);
            theDaughtersName2 = "";
            theDaughtersName3 = "";
            break;
         case 2:    
            theDaughtersName1 = theDecayChannel->GetDaughterName(0);
            theDaughtersName2 = theDecayChannel->GetDaughterName(1);
            theDaughtersName3 = "";
	    break;		
	 default:    
            theDaughtersName1 = theDecayChannel->GetDaughterName(0);
            theDaughtersName2 = theDecayChannel->GetDaughterName(1);
            theDaughtersName3 = theDecayChannel->GetDaughterName(2);
	    break;
	}

//	
//      Get the decay products List
//
     
     G4GeneralPhaseSpaceDecay thePhaseSpaceDecayChannel(theParentName,
                                                        theParentMass,
                                                        theBR,
                                                        theNumberOfDaughters,
                                                        theDaughtersName1,                  
		                                        theDaughtersName2,
		                                        theDaughtersName3);
     G4DecayProducts* theDecayProducts = thePhaseSpaceDecayChannel.DecayIt();
		                                        
//
//      Create the kinetic track List associated to the decay products
//
     G4LorentzRotation toMoving(Get4Momentum().boostVector());
     G4DynamicParticle* theDynamicParticle;
     G4double theFormationTime = 0.0;
     G4ThreeVector thePosition = this->GetPosition();
     G4LorentzVector momentum;
     G4KineticTrackVector* theDecayProductList = new G4KineticTrackVector;
     G4int dEntries = theDecayProducts->entries();
     for (G4int i=dEntries; i > 0; i--)
        {
         theDynamicParticle = theDecayProducts->PopProducts();
         theDefinition = theDynamicParticle->GetDefinition();
         momentum = toMoving*theDynamicParticle->Get4Momentum();
         theDecayProductList->push_back(new G4KineticTrack (theDefinition,
                                                         theFormationTime,
                                                         thePosition,
                                                         momentum));
         delete theDynamicParticle;
        }
     delete theDecayProducts;
     delete [] theCumActualWidth;
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
    sqrt(std::max((((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
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
    sqrt(std::max((((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
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
    sqrt(((mass * mass) - (G4KineticTrack_xmass1 + xmass) * (G4KineticTrack_xmass1 + xmass)) *
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








