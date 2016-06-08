// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KineticTrack.cc,v 1.3.8.1.2.1 1999/12/08 17:34:44 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
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
//   Default constructor
//

G4KineticTrack::G4KineticTrack() :
                theDefinition(NULL),
                theFormationTime(0.0),
                thePosition(0.0, 0.0, 0.0),
                the4Momentum(0.0,0.0,0.0,0.0),
                theInitialCoordinates(0.0,0.0,0.0,0.0),
                nChannels(0),
                theActualMass(0.0),            
                theActualWidth(NULL),            
                theDaughterMass(NULL),
                theDaughterWidth(NULL)
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
  theDaughterMass = NULL;
  theDaughterWidth = NULL;
 
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
 if (theDecayTable != NULL)
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

  theDaughterMass = NULL;
  theDaughterWidth = NULL;
  theActualWidth = NULL;
  G4bool * theDaughterIsShortLived = 0;
  
  if(nChannels!=0) theActualWidth = new G4double[nChannels];

  G4int index;
  for (index = nChannels - 1; index >= 0; index--)
     {
      G4VDecayChannel* theChannel = theDecayTable->GetDecayChannel(index);
      G4int nDaughters = theChannel->GetNumberOfDaughters();

      if (nDaughters == 2) 
         {
          G4double thePoleMass = theDefinition->GetPDGMass();
          G4double thePoleWidth = theChannel->GetBR();
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

              theActualMom = EvaluateCMMomentum(theActualMass, 
                                                theDaughterMass);   
              thePoleMom = EvaluateCMMomentum(thePoleMass, 
                                              theDaughterMass);

             }
          else if ( !theDaughterIsShortLived[0] && theDaughterIsShortLived[1] )   
             {

//              G4cout << G4endl << "Only the first of the " << nDaughters <<
//                              " decay products is stable!";

              theActualMom = IntegrateCMMomentum();
              thePoleMom = IntegrateCMMomentum(thePoleMass);
                
             }        
          else if ( theDaughterIsShortLived[0] && !theDaughterIsShortLived[1] )   
             {

//              G4cout << G4endl << "Only the second of the " << nDaughters <<
//                              " decay products is stable!";

//
//               Swap the content of the theDaughterMass and theDaughterWidth arrays!!!
//

              G4SwapObj(theDaughterMass, theDaughterMass + 1);
              G4SwapObj(theDaughterWidth, theDaughterWidth + 1);

              theActualMom = IntegrateCMMomentum();
              thePoleMom = IntegrateCMMomentum(thePoleMass);
                
             }        
          else if ( theDaughterIsShortLived[0] && theDaughterIsShortLived[1] )   
             {

//              G4cout << G4endl << "Both the " << nDaughters <<
//                              " decay products are resonances!";

              /* @@@@@@ code has to be implemented */

             }        

          G4double theMassRatio = thePoleMass / theActualMass;
          G4double theMomRatio = theActualMom / thePoleMom;
          G4int l = 1;   /* code has to be made more general */

//
//           Evaluate tha "actual widths"
//

          theActualWidth[index] = thePoleWidth * theMassRatio *
                                  pow(theMomRatio, (2 * l + 1)) *
                                  (1.2 / (0.2 + pow(theMomRatio, (2 * l))));

          delete [] theDaughterMass;
	  theDaughterMass = NULL;
          delete [] theDaughterWidth;
	  theDaughterWidth = NULL;
         }
      else
         {
          theActualWidth[index] = theChannel->GetBR();
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
 if (theActualWidth != NULL) delete [] theActualWidth;
 if (theDaughterMass != NULL) delete [] theDaughterMass;
 if (theDaughterWidth != NULL) delete [] theDaughterWidth;
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
     if (theActualWidth != NULL) delete [] theActualWidth;
     theActualWidth = new G4double[nChannels];
     for (i = 0; i < nChannels; i++) 
        {
         theActualWidth[i] = right.theActualWidth[i];
        }
     nChannels = right.GetnChannels();      
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
        }
     G4double r = theTotalActualWidth * G4UniformRand();
     G4ParticleDefinition* theDefinition = this->GetDefinition();
     G4DecayTable* theDecayTable = theDefinition->GetDecayTable();
     G4VDecayChannel* theDecayChannel;
     for (index = nChannels - 1; index >= 0; index--)
        {
         if (r < theCumActualWidth[index])
            {
             theDecayChannel = theDecayTable->GetDecayChannel(index);
             break; 
            }
        }
        
     G4String theParentName = theDecayChannel->GetParentName();
     G4double theParentMass = this->GetActualMass();
     G4double theBR = theActualWidth[index];
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
     return NULL;
    }
}



G4double G4KineticTrackIntegrandFunction1(G4double xmass)
{
 G4double result = 1.0;
// G4double mass = theActualMass;   /* the actual mass value */
// G4double mass1 = theDaughterMass[0];
// G4double mass2 = theDaughterMass[1];
// G4double gamma2 = theDaughtherWidth[1];
// result = (mass1 / (2 * mass)) *
//          sqrt(((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
//               ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))) *
//          twopi * (gamma2 / ((mass1 - xmass) * (mass2 - xmass) + ((gamma2 * gamma2) / 4)));    
 return result;
}

G4double G4KineticTrackIntegrandFunction2(G4double xmass)
{
 G4double result = 2.0;
// G4double mass = theDefinition->GetPDGMass();   /* the pole mass value */
// G4double mass1 = theDaughterMass[0];
// G4double mass2 = theDaughterMass[1];
// G4double gamma2 = theDaughtherWidth[1];
// result = (mass1 / (2 * mass)) *
//          sqrt(((mass * mass) - (mass1 + xmass) * (mass1 + xmass)) *
//               ((mass * mass) - (mass1 - xmass) * (mass1 - xmass))) *
//          twopi * (gamma2 / ((mass1 - xmass) * (mass2 - xmass) + ((gamma2 * gamma2) / 4)));    
//
 return result;
}



G4double G4KineticTrack::IntegrateCMMomentum() const
{
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = theActualMass - theDaughterMass[0];
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrackIntegrandFunction1);
 G4double theIntegralOverMass2 = theIntegrandMomentum.Simpson(theLowerLimit,
                                                              theUpperLimit,
                                                              nIterations);                                    

// G4cout << G4endl << "Inside IntegrateCMMomentum (actual mass case): ";
// G4cout << G4endl << "   Integration result = " << theIntegralOverMass2;

 return theIntegralOverMass2;
}

G4double G4KineticTrack::IntegrateCMMomentum(const G4double polemass) const
{
 G4double theLowerLimit = 0.0;
 G4double theUpperLimit = theActualMass - theDaughterMass[0];
 G4int nIterations = 100;
 G4SimpleIntegration theIntegrandMomentum(&G4KineticTrackIntegrandFunction2);
 G4double theIntegralOverMass2 = theIntegrandMomentum.Simpson(theLowerLimit,
                                                              theUpperLimit,
                                                              nIterations);                                    

// G4cout << G4endl << "Inside IntegrateCMMomentum (pole mass case): ";
// G4cout << G4endl << "   Integration result = " << theIntegralOverMass2;

 return theIntegralOverMass2;
}
