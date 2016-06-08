#ifndef G4QGSMSplitableHadron_h
#define G4QGSMSplitableHadron_h 1

#include "G4VSplitableHadron.hh"
#include "G4PartonVector.hh"
#include "G4MesonSplitter.hh"
#include "G4BaryonSplitter.hh"
#include "Randomize.hh"

// based on prototype by Maxim Komogorov
// Splitting into methods, and centralizing of model parameters HPW Feb 1999
// continued clean-up of interfaces and algorithms HPW 1999.
// Redesign of data structures and algorithms HPW Feb 1999

class G4QGSMSplitableHadron : public G4VSplitableHadron
{

  public:
      G4QGSMSplitableHadron();
      G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary);
      G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary, G4bool Direction);
      G4QGSMSplitableHadron(const G4Nucleon & aNucleon); 

      virtual ~G4QGSMSplitableHadron();
     
      const G4QGSMSplitableHadron & operator=(const G4QGSMSplitableHadron &right);

      virtual void SplitUp();
      virtual G4Parton * GetNextParton();
      virtual G4Parton * GetNextAntiParton();

  private:
      void InitParameters();
      void DiffractiveSplitUp();
      void SoftSplitUp();
      G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare);
      void GetValenceQuarkFlavors(const G4ParticleDefinition * aPart, 
                                  G4Parton *& Parton1, G4Parton *& Parton2);
      G4Parton * BuildSeaQuark(G4bool isAntiQuark, G4int aPDGCode, G4int nSeaPair);
      G4double SampleX(G4double anXmin, G4int nSea, G4int theTotalSea, G4double aBeta);
  
  private:
  // aggregated data
    G4bool Direction; // FALSE is target. - candidate for more detailed design. @@@@ HPW

    G4PartonVector Color;
    G4PartonVector AntiColor;   
  private:
    // associated classes
    G4MesonSplitter theMesonSplitter;
    G4BaryonSplitter theBaryonSplitter;
  
  private:
   // model parameters
    double alpha;
    double beta;
    double theMinPz; 
    double StrangeSuppress;
    double sigmaPt;
    double widthOfPtSquare; 
};

inline G4Parton* G4QGSMSplitableHadron::GetNextParton()
   {
   return Color.isEmpty()?0:Color.removeLast();
   }

inline G4Parton* G4QGSMSplitableHadron::GetNextAntiParton()
   {
   return AntiColor.isEmpty()?0:AntiColor.removeFirst();
   }

#endif


