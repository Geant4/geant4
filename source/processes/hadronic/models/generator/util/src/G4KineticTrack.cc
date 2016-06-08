// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KineticTrack.cc,v 1.10 2000/12/07 14:38:21 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// $Id: G4KineticTrack.cc,v 1.0 1998/05/20
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
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

//
// Some static clobal for integration
//

static G4double  G4KineticTrack_Gmass, G4KineticTrack_Gmass1, G4KineticTrack_Gmass2,G4KineticTrack_Ggamma1,G4KineticTrack_Ggamma2,G4KineticTrack_xmass1;

//
//   Default constructor
//

G4KineticTrack::G4KineticTrack() :
                theDefinition(0),
                theFormationTime(0.0),
                thePosition(0.0, 0.0, 0.0),
                the4Momentum(0.0,0.0,0.0,0.0),
                theInitialCoordinates(0.0,0.0,0.0,0.0),
                nChannels(0),
                theActualMass(0.0),            
                theActualWidth(0),            
                theDaughterMass(0),
                theDaughterWidth(0)
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

G4KineticTrack::G4KineticTrack(const G4KineticTrack &right)
{
 G4int i;
 theDefinition = right.GetDefinition();
 theFormationTime = right.GetFormationTime();
 thePosition = right.GetPosition();
 the4Momentum = right.Get4Momentum();
 theInitialCoordinates = right.GetInitialCoordinates();
 nChannels = right.GetnChannels();
 theActualMass = right.GetActualMass();
 theActualWidth = new G4double[nChannels];
 for (i = 0; i < nChannels; i++)
  {
    theActualWidth[i] = right.theActualWidth[i];
  }
  theDaughterMass = 0;
  theDaughterWidth = 0;
 
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
                the4Momentum(a4Momentum)
{
//
//   Initialize the InitialCoordinates 4vector
//

 theInitialCoordinates.setVect(thePosition);
 theInitialCoordinates.setT(theFormationTime); 

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
      if (nDaughters == 2) 
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
               G4KineticTrack_Gmass = theActualMass;
               G4KineticTrack_Gmass1 = theDaughterMass[0];
               G4KineticTrack_Gmass2 = theDaughterMass[1];
              G4KineticTrack_Ggamma2 = theDaughterWidth[1];
              theActualMom = IntegrateCMMomentum();

               G4KineticTrack_Gmass = thePoleMass;
              thePoleMom = IntegrateCMMomentum(thePoleMass);
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
               G4KineticTrack_Gmass = theActualMass;
               G4KineticTrack_Gmass1 = theDaughterMass[0];
               G4KineticTrack_Gmass2 = theDaughterMass[1];
              G4KineticTrack_Ggamma2 = theDaughterWidth[1];
              theActualMom = IntegrateCMMomentum();

               G4KineticTrack_Gmass = thePoleMass;
              thePoleMom = IntegrateCMMomentum(thePoleMass);
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
               G4KineticTrack_Gmass1 = theDaughterMass[0];
              G4KineticTrack_Ggamma1 = theDaughterWidth[0];
               G4KineticTrack_Gmass2 = theDaughterMass[1];
              G4KineticTrack_Ggamma2 = theDaughterWidth[1];
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

          G4double theMassRatio = thePoleMass / theActualMass;
          G4double theMomRatio = theActualMom / thePoleMom;
          G4int l = 0;   /* code has to be made more general */

//
//           Evaluate tha "actual widths"
//

          theActualWidth[index] = thePoleWidth * theMassRatio *
                                  pow(theMomRatio, (2 * l + 1)) *
                                  (1.2 / (1+ 0.2*pow(theMomRatio, (2 * l))));

	  //	  	  cout << " LB The Mass Ratio   = " << theMassRatio << G4endl;
	  //	  	  cout << " LB The Mom. Ratio   = " << theMomRatio << G4endl;
	  //	  	  cout << " LB The Pole Width   = " << thePoleWidth << G4endl;
	  //	  	  cout << " LB The Actual Width = " << theActualWidth[index] << "index: " << index << G4endl;

          delete [] theDaughterMass;
	  theDaughterMass = 0;
          delete [] theDaughterWidth;
	  theDaughterWidth = 0;
	  delete [] theDaughterIsShortLived;
          theDaughterIsShortLived = 0;
         }
      else
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
     thePosition = right.GetPosition();
     the4Momentum = right.Get4Momentum();  
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
     G4VDecayChannel* theDecayChannel;
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
     G4LorentzVector the4Momentum;
     G4KineticTrackVector* theDecayProductList = new G4KineticTrackVector;
     G4int dEntries = theDecayProducts->entries();
     for (G4int i=dEntries; i > 0; i--)
        {
         theDynamicParticle = theDecayProducts->PopProducts();
         theDefinition = theDynamicParticle->GetDefinition();
         the4Momentum = toMoving*theDynamicParticle->Get4Momentum();
         theDecayProductList->insert(new G4KineticTrack (theDefinition,
                                                         theFormationTime,
                                                         thePosition,
                                                         the4Momentum));
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

G4double G4KineticTrack_IntegrandFunction1(G4double xmass)
{
  // G4double mass = theActualMass;   /* the actual mass value */
 G4double mass =  G4KineticTrack_Gmass;   /* the actual mass value */
 // G4double mass1 = theDaughterMass[0];
 G4double mass1 =  G4KineticTrack_Gmass1;
 // G4double mass2 = theDaughterMass[1];
 G4double mass2 =  G4KineticTrack_Gmass2;
 // G4double gamma2 = theDaughterWidth[1];
 G4double gamma2 = G4KineticTrack_Ggamma2;
 G4KineticTrack wei;
 G4double result = (1. / (2 * mass)) *
          sqrt(((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
               ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))) *
          wei.BrWig(gamma2, mass2, xmass);
 return result;
}

G4double G4KineticTrack_IntegrandFunction2(G4double xmass)
{
  // G4double mass = theDefinition->GetPDGMass();   /* the pole mass value */
 G4double mass =  G4KineticTrack_Gmass;   /* the actual mass value */
  // G4double mass1 = theDaughterMass[0];
 G4double mass1 =  G4KineticTrack_Gmass1;
  // G4double mass2 = theDaughterMass[1];
 G4double mass2 =  G4KineticTrack_Gmass2;
  // G4double gamma2 = theDaughterWidth[1];
 G4double gamma2 = G4KineticTrack_Ggamma2;
 G4KineticTrack wei;
 G4double result = (1. / (2 * mass)) *
          sqrt(((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
               ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))) *
          wei.BrWig(gamma2, mass2, xmass);
 return result;
}

G4double G4KineticTrack_IntegrandFunction3(G4double xmass)
{
 G4double mass =  G4KineticTrack_Gmass;   /* the actual mass value */
 G4double mass1 =  G4KineticTrack_Gmass1;
 G4double gamma1 = G4KineticTrack_Ggamma1;
 G4double mass2 =  G4KineticTrack_Gmass2;
 G4double gamma2 = G4KineticTrack_Ggamma2;
 G4KineticTrack wei;
 G4double result = (1. / (2 * mass)) *
          sqrt(((mass * mass) - (G4KineticTrack_xmass1 + xmass) * (G4KineticTrack_xmass1 + xmass)) *
               ((mass * mass) - (G4KineticTrack_xmass1 - xmass) * (G4KineticTrack_xmass1 - xmass))) *
          wei.BrWig(gamma2, mass2, xmass);
 return result;
}

G4double G4KineticTrack_IntegrandFunction4(G4double xmass)
{
 G4double mass =  G4KineticTrack_Gmass;
 G4double mass1 =  G4KineticTrack_Gmass1;
 G4double gamma1 = G4KineticTrack_Ggamma1;
 G4double mass2 =  G4KineticTrack_Gmass2;
 G4double gamma2 = G4KineticTrack_Ggamma2;
 G4KineticTrack wei;
 G4double G4KineticTrack_xmass1=xmass;
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = mass - xmass;
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrack_IntegrandFunction3);
 G4double result = wei.BrWig(gamma1, mass1, xmass)*
                   theIntegrandMomentum.Simpson(theLowerLimit,
						theUpperLimit,
						nIterations);
 return result;
}

G4double G4KineticTrack::IntegrateCMMomentum() const
{
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = theActualMass - theDaughterMass[0];
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrack_IntegrandFunction1);
 G4double theIntegralOverMass2;
 if(theUpperLimit>0.)
  theIntegralOverMass2 = theIntegrandMomentum.Simpson(theLowerLimit,
                                                              theUpperLimit,
                                                              nIterations);                                    
 else
   theIntegralOverMass2 = 0.;

 // G4cout << G4endl << "Inside IntegrateCMMomentum (actual mass case): ";
 // G4cout << G4endl << "   Integration result = " << theIntegralOverMass2 << G4endl ;
 return theIntegralOverMass2;
}

G4double G4KineticTrack::IntegrateCMMomentum(const G4double polemass) const
{
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = polemass - theDaughterMass[0];
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrack_IntegrandFunction1);
 G4double theIntegralOverMass2;
 if(theUpperLimit>0.)
  theIntegralOverMass2 = theIntegrandMomentum.Simpson(theLowerLimit,
                                                              theUpperLimit,
                                                              nIterations);  
 else
   theIntegralOverMass2 = 0.;

 // G4cout << G4endl << "Inside IntegrateCMMomentum (pole mass case): ";
 // G4cout << G4endl << "   Integration result = " << theIntegralOverMass2 << G4endl;
 // G4cout << G4endl << " Lower / upper limit " << theLowerLimit << " / " << theUpperLimit << G4endl;

 return theIntegralOverMass2;
}

G4double G4KineticTrack::IntegrateCMMomentum2() const
{
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = theActualMass;
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrack_IntegrandFunction4);
 G4double theIntegralOverMass2;
 if(theUpperLimit>0.)
  theIntegralOverMass2 = theIntegrandMomentum.Simpson(theLowerLimit,
                                                              theUpperLimit,
                                                              nIterations);  
 else
   theIntegralOverMass2 = 0.;

 // G4cout << G4endl << "Inside IntegrateCMMomentum 2: ";
 // G4cout << G4endl << "   Integration result = " << theIntegralOverMass2 << G4endl;

 return theIntegralOverMass2;
}








