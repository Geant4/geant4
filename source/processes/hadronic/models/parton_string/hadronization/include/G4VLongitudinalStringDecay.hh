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
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
// -----------------------------------------------------------------------------
#ifndef G4VLongitudinalStringDecay_h
#define G4VLongitudinalStringDecay_h 1

#include "G4HadronicInteraction.hh"
#include "G4VStringFragmentation.hh"
#include "G4DynamicParticle.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4HadronBuilder.hh"
#include <vector>

//*****************************************************************************

class G4FragmentingString;

class G4VLongitudinalStringDecay : public G4HadronicInteraction
{
  public:

    G4VLongitudinalStringDecay(const G4String& name = "StringDecay");
    virtual ~G4VLongitudinalStringDecay();

    G4HadFinalState *ApplyYourself(const G4HadProjectile&, G4Nucleus&) final;

  private:
    // not implemented to protect/forbid use
    G4VLongitudinalStringDecay(const G4VLongitudinalStringDecay &right);
    const G4VLongitudinalStringDecay & operator=(const G4VLongitudinalStringDecay &right);
    G4bool operator==(const G4VLongitudinalStringDecay &right) const;
    G4bool operator!=(const G4VLongitudinalStringDecay &right) const;

  public:
    virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString)=0;

    void AddNewParticles();
    void EraseNewParticles();
    //struct DeleteString { void operator()(G4ExcitedString* aS){delete aS;} };

    // To set minimal mass of a string. The string with mass above the minimal mass can fragment.
    void SetMinMasses();
    void SetMinimalStringMass(const G4FragmentingString  * const string);
    void SetMinimalStringMass2(const G4double aValue);

  protected: 
    // For changing Mass Cut used for selection of very small mass strings
    virtual void SetMassCut(G4double aValue);
    G4double GetMassCut();

    // For handling a string with very low mass
    G4KineticTrackVector * ProduceOneHadron(const G4ExcitedString * const theString);

    // To store created quarks or 2 last hadrons
    typedef std::pair<G4ParticleDefinition*, G4ParticleDefinition*> pDefPair;

    // For creation of hadrons from given quark pair 
    typedef G4ParticleDefinition * (G4HadronBuilder::*Pcreate)
		   		     (G4ParticleDefinition*, G4ParticleDefinition*);

    // Used by ProduceOneHadron method for estimation of lowest possible mass of
    // given quark system -- string.
    G4double PossibleHadronMass(const G4FragmentingString * const string,
	 	               Pcreate build=0, pDefPair * pdefs=0);

    G4ParticleDefinition* FindParticle(G4int Encoding); 

    // For decision on continue or stop string fragmentation
    virtual G4bool StopFragmenting(const G4FragmentingString  * const string)=0;
    virtual G4bool IsItFragmentable(const G4FragmentingString * const string)=0;

    // If a string can not fragment, make last break into 2 hadrons
    virtual G4bool SplitLast(G4FragmentingString * string, 
	   	             G4KineticTrackVector * LeftVector,
		             G4KineticTrackVector * RightVector)=0;

    virtual void Sample4Momentum(G4LorentzVector* Mom,     G4double Mass, 
                                 G4LorentzVector* AntiMom, G4double AntiMass, 
                                 G4double InitialMass)=0; 

    // If a string can fragment, do the following:

    // Make a copy of a string
    G4ExcitedString *CopyExcited(const G4ExcitedString& string);

    // Produce a hadron at Splitup of the string
    virtual G4KineticTrack * Splitup(G4FragmentingString *string,
                                     G4FragmentingString *&newString)=0;

    // The hadron can be producet at QuarkSplitup or DiQuarkSplitup
    virtual G4ParticleDefinition * QuarkSplitup(G4ParticleDefinition* decay,
	        	   		        G4ParticleDefinition *&created);

    virtual G4ParticleDefinition * DiQuarkSplitup(G4ParticleDefinition* decay,
		   			          G4ParticleDefinition *&created)=0;

    // All of them are going through quak-antiquark pair creation					
    pDefPair CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks=true);

  public:
    // For a pair it is needed:
    G4int SampleQuarkFlavor(void);
    G4ThreeVector SampleQuarkPt(G4double ptMax=-1.); // -1. no limit on maxpt.

  protected:
    // For determination of kinematical properties of the created hadron
    virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron,      
                                         G4FragmentingString * string,        
                                         G4FragmentingString * newString  )=0;

    virtual G4double GetLightConeZ(G4double zmin, G4double zmax, 
                                   G4int PartonEncoding,  
                                   G4ParticleDefinition* pHadron, 
                                   G4double Px, G4double Py       ) = 0;      

    void CalculateHadronTimePosition(G4double theInitialStringMass, 
                                     G4KineticTrackVector *);

    // Used for some test purposes 
    void ConstructParticle();

    G4ParticleDefinition* CreateHadron(G4int id1, G4int id2, 
                                       G4bool theGivenSpin, G4int theSpin); 

  public:
    void SetSigmaTransverseMomentum(G4double aQT);
    void SetStrangenessSuppression(G4double aValue);
    void SetDiquarkSuppression(G4double aValue);
    void SetDiquarkBreakProbability(G4double aValue);

    void SetSpinThreeHalfBarionProbability(G4double aValue);
   
    void SetScalarMesonMixings( std::vector<G4double> aVector);
    void SetVectorMesonMixings( std::vector<G4double> aVector);

    void SetStringTensionParameter(G4double aValue);            

    void SetProbCCbar(G4double aValue);
    void SetProbEta_c(G4double aValue);
    void SetProbBBbar(G4double aValue);
    void SetProbEta_b(G4double aValue);

  protected:  
    G4double GetDiquarkSuppress()	{return DiquarkSuppress;};
    G4double GetDiquarkBreakProb()	{return DiquarkBreakProb;};
    G4double GetStrangeSuppress()	{return StrangeSuppress;};
    G4int    GetClusterLoopInterrupt()   {return ClusterLoopInterrupt;};

    G4double GetProbCCbar(){return ProbCCbar;};
    G4double GetProbEta_c(){return ProbEta_c;};
    G4double GetProbBBbar(){return ProbBBbar;};
    G4double GetProbEta_b(){return ProbEta_b;};

    G4double GetStringTensionParameter() {return Kappa;};

  protected:  
    G4double  MassCut;
    G4double  SigmaQT;          // sigma_q_t of quark/hadron transverse momentum distribution parameter 
    G4double  DiquarkSuppress;  // Diquark suppression parameter  
    G4double  DiquarkBreakProb; // Diquark breaking probability, qq->h+qq'
    G4double  StrangeSuppress ;
    G4int     StringLoopInterrupt;
    G4int     ClusterLoopInterrupt;

    G4HadronBuilder *hadronizer;

    std::vector<G4double> pspin_meson;
    G4double pspin_barion;
    std::vector<G4double> vectorMesonMix;
    std::vector<G4double> scalarMesonMix;
   
    G4double ProbCCbar;  // Probability of C-Cbar pair creation
    G4double ProbEta_c;  // Mixing of Eta_c and J/Psi

    G4double ProbBBbar;  // Probability of B-Bbar pair creation
    G4double ProbEta_b;  // Mixing of Eta_b and Ipsilon_b

    G4double ProbCB;     // = ProbCCbar + ProbBBbar

    G4double MaxMass;

    G4bool   PastInitPhase;

    G4double Kappa;      // String tension parameter

    std::vector<G4ParticleDefinition *> NewParticles;

  public:
    // ------ For estimation of a minimal string mass ---------------
    G4double Mass_of_light_quark;
    G4double Mass_of_s_quark;
    G4double Mass_of_c_quark;
    G4double Mass_of_b_quark;
    G4double Mass_of_string_junction;

    G4double minMassQQbarStr[5][5];
    G4double minMassQDiQStr[5][5][5];

    // ------ An estimated minimal string mass ----------------------
    G4double MinimalStringMass;
    G4double MinimalStringMass2;

    G4int Qcharge[5];  // quark charges
    G4int          Meson[5][5][7];
    G4double MesonWeight[5][5][7];

    G4int          Baryon[5][5][5][4];
    G4double BaryonWeight[5][5][5][4];

    G4double Prob_QQbar[5];

    G4int DecayQuark;
    G4int NewQuark;
    /*
    G4double FFq2q[5][5][2];
    G4double FFq2qq[5][15][2];
    G4double FFqq2q[15][5][2];
    G4double FFqq2qq[15][5][2];
    */

    // ------ To improve the code structure
    G4ParticleDefinition * FS_LeftHadron[350], * FS_RightHadron[350];
    G4double FS_Weight[350];
    G4int NumberOf_FS;
};

//******************************************************************************
// Class G4VLongitudinalStringDecay 

#endif

