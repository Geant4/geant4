// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VLongitudinalStringDecay.hh,v 1.4 1998/12/01 15:40:47 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
// Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
// -----------------------------------------------------------------------------
#ifndef G4VLongitudinalStringDecay_h
#define G4VLongitudinalStringDecay_h 1
#include "G4VStringFragmentation.hh"
#include "G4DynamicParticle.hh"

//**********************************************************************************************

class G4VLongitudinalStringDecay: public G4VStringFragmentation 
   {
public:
   G4VLongitudinalStringDecay();     
  ~G4VLongitudinalStringDecay();

private:
//  G4VLongitudinalStringDecay(const G4VLongitudinalStringDecay &right);
//  const G4VLongitudinalStringDecay & operator=(const G4VLongitudinalStringDecay &right);

   int operator==(const G4VLongitudinalStringDecay &right) const;
   int operator!=(const G4VLongitudinalStringDecay &right) const;

public:
   G4KineticTrackVector* FragmentString(const G4ExcitedString& theString);
   G4bool FragmentString(G4KineticTrackVector* aHadrons, const G4ExcitedString* theString);
   G4KineticTrackVector* DecayResonans (G4KineticTrackVector* aHadrons);

   G4int SampleQuarkFlavor(void);
   void  SampleQuarkPt(G4double* thePx, G4double* thePy);

protected:
   // Additional protected declarations 
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py) = 0;      
   G4ParticleDefinition* CreateHadron(G4int id1, G4int id2, G4bool theGivenSpin, G4int theSpin); 
   void CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector *);
   G4ParticleDefinition* FindParticle(G4int Encoding); 
   void Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass); 

   // Additional Implementation Declarations

private:  
   G4double  MassCut;
   G4double  ClusterMass;
   G4double  SigmaQT;          // sigma_q_t is quark transverse momentum distribution parameter 
   G4double  DiquarkSuppress;  // is Diquark suppression parameter  
   G4double  DiquarkBreakProb; // is Diquark breaking probability 
   G4double  SmoothParam;      // model parameter
   G4double  StrangeSuppress ;
   G4int     StringLoopInterrupt;
   G4int     ClusterLoopInterrupt;

   void ConstructParticle();
};

//**********************************************************************************************
// Class G4VLongitudinalStringDecay 
#endif


