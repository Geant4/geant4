// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VLongitudinalStringDecay.hh,v 1.3 1999/05/20 15:22:31 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
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
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

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
   G4double GetDiquarkSuppress()	{return DiquarkSuppress;};
   G4double GetDiquarkBreakProb()	{return DiquarkBreakProb;};
   G4double GetStrangeSuppress()	{return StrangeSuppress;};
   G4double GetClusterMass()		{return ClusterMass;};
   G4int    GetClusterLoopInterrupt()   {return ClusterLoopInterrupt;};
   
   G4ParticleDefinition* CreateHadron(G4int id1, G4int id2, G4bool theGivenSpin, G4int theSpin); 
   void Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass); 

protected:
   // Additional protected declarations 
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py) = 0;      
   void CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector *);
   G4ParticleDefinition* FindParticle(G4int Encoding); 

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
   
   
   
class SimpleString
{
private:
	class SideOfString;

public:
	SimpleString(const G4ExcitedString& excitedString, G4VLongitudinalStringDecay * stringdecay)
	: left(new SideOfString(excitedString.GetLeftParton()->GetPDGcode(),excitedString.Get4Momentum().mag())),
	  right(new SideOfString(excitedString.GetRightParton()->GetPDGcode(),excitedString.Get4Momentum().mag())),
	  MassSquare(excitedString.Get4Momentum().mag2()),
	  theStringDecay(stringdecay),
	  decay(0), stable(0), Side(0)
	{};
	~SimpleString()
	{
		delete right;
		delete left;
	};
	SideOfString * Left(){return left;};
	SideOfString * Right(){return right;};
	SideOfString * Stable(){return stable;};
	SideOfString * Decay(){return decay;};
	const G4double MassSquared(){return MassSquare;};
	G4ParticleDefinition * Splitup(G4int & NewDecayEncoding);
	G4bool SplitLast(G4KineticTrackVector * LeftVector,G4KineticTrackVector * RightVector);
	G4int GetDecayDirection() { return Side;};
	void update();

private:
	SimpleString(){cout << " error calling bad ctor for class SimpleString"<< endl;};
	
	SideOfString * left, * right;
	SideOfString * decay, * stable;
	G4double MassSquare;
	G4VLongitudinalStringDecay *theStringDecay;
	G4int Side;
	class SideOfString
	{      
	public:   
	   G4int    Encoding(){return theEncoding;};
	   G4double w(){return thew;};
	   G4double Px(){return thepx;};
	   G4double Py(){return thepy;};
	   void decreasew(G4double deltaw){thew -= deltaw;};
	   void setEncoding(G4int aEncoding){theEncoding=aEncoding;};
	   void setpxpy(G4double px,G4double py){thepx=px; thepy=py;};

	   SideOfString(G4int aEncoding, G4double aw)
	   { 
   		theEncoding = aEncoding; 
   		thew = aw; 
   		thepx = thepy = 0; 
	   }
	   ;
	   SideOfString(void) {};
	//   void Init(G4int aEncoding, G4double aw)   { 
	//    	theEncoding = aEncoding; 
	//    	thew = aw; 
	//    	thepx = thepy = 0; 
	//    };

	private:
	   G4int theEncoding;
	   G4double thew, thepx, thepy;

	};
};

   G4KineticTrackVector * LightFragmentationTest(const G4ExcitedString& theString);
   G4bool StopFragmenting(SimpleString& string);
   G4double IsFragmentable(SimpleString & theString);
   G4double MinFragmentationMass(SimpleString & theString,
				G4ParticleDefinition*& Hadron1,
				G4ParticleDefinition*& Hadron2);
   G4KineticTrack * SplitEandP(G4ParticleDefinition * pHadron,
                               SimpleString * string, G4int newDecayEncoding);

};


//**********************************************************************************************
// Class G4VLongitudinalStringDecay 
#endif


