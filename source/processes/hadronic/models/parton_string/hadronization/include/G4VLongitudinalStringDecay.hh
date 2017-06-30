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
// $Id: G4VLongitudinalStringDecay.hh 102717 2017-02-20 10:37:13Z gcosmo $
// Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
// -----------------------------------------------------------------------------
#ifndef G4VLongitudinalStringDecay_h
#define G4VLongitudinalStringDecay_h 1
#include "G4VStringFragmentation.hh"
#include "G4DynamicParticle.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4HadronBuilder.hh"

class G4FragmentingString;
//**************************************************************************************

class G4VLongitudinalStringDecay 
   {
public:
            G4VLongitudinalStringDecay();     
   virtual ~G4VLongitudinalStringDecay();

private:
   // not implemented to protect/forbid use
   G4VLongitudinalStringDecay(const G4VLongitudinalStringDecay &right);
   const G4VLongitudinalStringDecay & operator=(const G4VLongitudinalStringDecay &right);
   int operator==(const G4VLongitudinalStringDecay &right) const;
   int operator!=(const G4VLongitudinalStringDecay &right) const;

public:
   virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString)=0;

protected: 

// For changing Mass Cut used for selection of very small mass strings
   virtual void SetMassCut(G4double aValue);
   G4double GetMassCut();

// For handling a string with very low mass
   G4KineticTrackVector * LightFragmentationTest(const G4ExcitedString * const theString);

// To store created quarks or 2 last hadrons
   typedef std::pair<G4ParticleDefinition*, G4ParticleDefinition*> pDefPair;

// For creation of hadrons from given quark pair 
   typedef G4ParticleDefinition * (G4HadronBuilder::*Pcreate)
		   			(G4ParticleDefinition*, G4ParticleDefinition*);

//-----------------------------------------------------------------------------
// Used by LightFragmentationTest for estimation of lowest possible mass of
// given quark system
   G4double FragmentationMass(const G4FragmentingString * const string,
		              Pcreate build=0,
		              pDefPair * pdefs=0);

   G4ParticleDefinition* FindParticle(G4int Encoding); 

   virtual void Sample4Momentum(G4LorentzVector* Mom,     G4double Mass, 
                                G4LorentzVector* AntiMom, G4double AntiMass, 
                                G4double InitialMass)=0; 
//-----------------------------------------------------------------------------
// For decision on continue or stop string fragmentation
   virtual G4bool StopFragmenting(const G4FragmentingString  * const string)=0;
   virtual G4bool IsFragmentable(const G4FragmentingString * const string)=0;

// If a string can not fragment, make last break into 2 hadrons
   virtual G4bool SplitLast(G4FragmentingString * string, 
		    G4KineticTrackVector * LeftVector,
		    G4KineticTrackVector * RightVector)=0;
//-----------------------------------------------------------------------------

// If a string fragments, do the following

// To make a copy of a string
   G4ExcitedString *CopyExcited(const G4ExcitedString& string);

   G4ParticleDefinition * QuarkSplitup(G4ParticleDefinition* decay,
		   		       G4ParticleDefinition *&created);

   virtual G4ParticleDefinition * DiQuarkSplitup(G4ParticleDefinition* decay,
		   			 G4ParticleDefinition *&created)=0;
					
   pDefPair CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks=true);

public:
//   used by G4VKinkyStringDecy..
   G4int SampleQuarkFlavor(void);
   G4ThreeVector SampleQuarkPt(G4double ptMax=-1.); // -1. no limit on maxpt.

protected:

//-----------------------------------------------------------------------------
// For determination of kinematical properties of created hadron
//   virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron,    
//                                        G4FragmentingString * string  )=0; 
   virtual G4KineticTrack * Splitup(G4FragmentingString *string,              // Uzhi 28 June 2016
                            G4FragmentingString *&newString)=0;

   virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron,      
                                        G4FragmentingString * string,        
                                        G4FragmentingString * newString  )=0;

   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, 
                                  G4int PartonEncoding,  
                                  G4ParticleDefinition* pHadron, 
                                  G4double Px, G4double Py       ) = 0;      

   void CalculateHadronTimePosition(G4double theInitialStringMass, 
                                    G4KineticTrackVector *);

// Used for some test purposes ------------------------------------------------ 
   void ConstructParticle();

   G4ParticleDefinition* CreateHadron(G4int id1, G4int id2, 
                                      G4bool theGivenSpin, G4int theSpin); 
   
//-----------------------------------------------------------------------------
public:

   G4KineticTrackVector* DecayResonans (G4KineticTrackVector* aHadrons);

   void SetSigmaTransverseMomentum(G4double aQT);
   void SetStrangenessSuppression(G4double aValue);
   void SetDiquarkSuppression(G4double aValue);
   void SetDiquarkBreakProbability(G4double aValue);

   void SetVectorMesonProbability(G4double aValue);
   void SetSpinThreeHalfBarionProbability(G4double aValue);
   
   void SetScalarMesonMixings( std::vector<G4double> aVector);
   void SetVectorMesonMixings( std::vector<G4double> aVector);

   void SetStringTensionParameter(G4double aValue);            

//   G4LorentzVector RestMomentum;                 // Uzhi June 2016
//private:
protected:  
   G4double GetDiquarkSuppress()	{return DiquarkSuppress;};
   G4double GetDiquarkBreakProb()	{return DiquarkBreakProb;};
   G4double GetStrangeSuppress()	{return StrangeSuppress;};
   G4double GetClusterMass()		{return ClusterMass;};
   G4int    GetClusterLoopInterrupt()   {return ClusterLoopInterrupt;};

   G4double GetStringTensionParameter() {return Kappa;};
   
//private:
protected:  
   G4double  MassCut;
   G4double  ClusterMass;
   G4double  SigmaQT;          // sigma_q_t is quark transverse momentum distribution parameter 
   G4double  DiquarkSuppress;  // is Diquark suppression parameter  
   G4double  DiquarkBreakProb; // is Diquark breaking probability 
   G4double  SmoothParam;      // model parameter
   G4double  StrangeSuppress ;
   G4int     StringLoopInterrupt;
   G4int     ClusterLoopInterrupt;
   G4HadronBuilder *hadronizer;

   G4double pspin_meson;
   G4double pspin_barion;
   std::vector<G4double> vectorMesonMix;
   std::vector<G4double> scalarMesonMix;
   
   G4bool    PastInitPhase;

   G4double Kappa; // String tension parameter

//   G4double MinFragmentationMass(G4ExcitedString * theString,
//				G4ParticleDefinition*& Hadron1,
//				G4ParticleDefinition*& Hadron2);

};

//*************************************************************************************
// Class G4VLongitudinalStringDecay 
#endif
