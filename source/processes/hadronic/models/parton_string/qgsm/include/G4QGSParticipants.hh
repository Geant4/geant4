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
#ifndef G4QGSParticipants_h
#define G4QGSParticipants_h 1

#include "Randomize.hh"
#include "G4VParticipants.hh"
#include "G4Nucleon.hh"
#include "G4InteractionContent.hh"
#include "G4QGSDiffractiveExcitation.hh"
#include "G4SingleDiffractiveExcitation.hh"
#include "G4PartonPair.hh" 
#include "G4QGSMSplitableHadron.hh" 
#include "G4V3DNucleus.hh"

#include "G4VSplitableHadron.hh"                 // Uzhi

#include "G4Reggeons.hh"
#include "G4QuarkExchange.hh"

class G4QGSParticipants : public G4VParticipants
{
public:
	G4QGSParticipants();
	G4QGSParticipants(const G4QGSParticipants &right);
	const G4QGSParticipants & operator=(const G4QGSParticipants &right);
	virtual ~G4QGSParticipants();

	int operator==(const G4QGSParticipants &right) const;
	int operator!=(const G4QGSParticipants &right) const;

	virtual void DoLorentzBoost(G4ThreeVector aBoost)
	{
                theCurrentVelocity = -aBoost;
		if(theNucleus) theNucleus->DoLorentzBoost(aBoost);
		theBoost = aBoost;
	}

	G4PartonPair* GetNextPartonPair();
	void BuildInteractions(const G4ReactionProduct  &thePrimary);
	void StartPartonPairLoop();

private:
    G4V3DNucleus* GetTargetNucleus() const;
    G4V3DNucleus* GetProjectileNucleus() const;

    void PrepareInitialState( const G4ReactionProduct& thePrimary );
    void GetList( const G4ReactionProduct& thePrimary );

    void StoreInvolvedNucleon();              
    void ReggeonCascade();
    G4bool PutOnMassShell();
    void GetResiduals();                  

    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;

    G4bool ComputeNucleusProperties( G4V3DNucleus* nucleus, G4LorentzVector& nucleusMomentum, 
                                     G4LorentzVector& residualMomentum, G4double& sumMasses,   
                                     G4double& residualExcitationEnergy, G4double& residualMass,
                                     G4int& residualMassNumber, G4int& residualCharge );
    // Utility methods used by PutOnMassShell.

    G4bool GenerateDeltaIsobar( const G4double sqrtS, const G4int numberOfInvolvedNucleons,
                                G4Nucleon* involvedNucleons[], G4double& sumMasses );

    G4bool SamplingNucleonKinematics( G4double averagePt2, const G4double maxPt2,
                                      G4double dCor, G4V3DNucleus* nucleus, 
                                      const G4LorentzVector& pResidual, 
                                      const G4double residualMass, const G4int residualMassNumber,
                                      const G4int numberOfInvolvedNucleons,
                                      G4Nucleon* involvedNucleons[], G4double& mass2 );

    G4bool CheckKinematics( const G4double sValue, const G4double sqrtS, 
                            const G4double projectileMass2, const G4double targetMass2,
                            const G4double nucleusY, const G4bool isProjectileNucleus,
                            const G4int numberOfInvolvedNucleons, G4Nucleon* involvedNucleons[],
                            G4double& targetWminus, G4double& projectileWplus, G4bool& success );

    G4bool FinalizeKinematics( const G4double w, const G4bool isProjectileNucleus, 
                               const G4LorentzRotation& boostFromCmsToLab,
                               const G4double residualMass, const G4int residualMassNumber,
                               const G4int numberOfInvolvedNucleons, 
                               G4Nucleon* involvedNucleons[],
	                       G4LorentzVector& residual4Momentum );

    void CreateStrings();

private:
    // Set parameters of nuclear destruction
    void SetCofNuclearDestruction( const G4double aValue );
    void SetR2ofNuclearDestruction( const G4double aValue );

    void SetExcitationEnergyPerWoundedNucleon( const G4double aValue );

    void SetDofNuclearDestruction( const G4double aValue );
    void SetPt2ofNuclearDestruction( const G4double aValue );
    void SetMaxPt2ofNuclearDestruction( const G4double aValue );

    // Get parameters of nuclear destruction
    G4double GetCofNuclearDestruction();
    G4double GetR2ofNuclearDestruction();

    G4double GetExcitationEnergyPerWoundedNucleon();

    G4double GetDofNuclearDestruction();
    G4double GetPt2ofNuclearDestruction();
    G4double GetMaxPt2ofNuclearDestruction();

protected:
	virtual G4VSplitableHadron* SelectInteractions(const G4ReactionProduct  &thePrimary);
	void SplitHadrons(); 
	void PerformSoftCollisions();
	void PerformDiffractiveCollisions();
        G4bool DeterminePartonMomenta();

protected:
	struct DeleteInteractionContent {void operator()(G4InteractionContent*aC){delete aC;}};
	std::vector<G4InteractionContent*> theInteractions;
	struct DeleteSplitableHadron{void operator()(G4VSplitableHadron*aS){delete aS;}};
	std::vector<G4VSplitableHadron*>   theTargets;
	struct DeletePartonPair{void operator()(G4PartonPair*aP){delete aP;}};
	std::vector<G4PartonPair*>   thePartonPairs;

	G4QuarkExchange 	      theQuarkExchange;         // Uzhi 20 Oct. 2016
	G4SingleDiffractiveExcitation theSingleDiffExcitation;
	G4QGSDiffractiveExcitation    theDiffExcitaton;
	G4int ModelMode;

	G4ThreeVector theBoost;
	G4double SampleX(G4double anXmin, G4int nSea, G4int theTotalSea, G4double aBeta);
protected:
	// model parameters HPW
	enum  { SOFT, DIFFRACTIVE };
	enum  { ALL, WITHOUT_R, NON_DIFF };   // Interaction modes
	enum  { PrD, TrD, DD, NonD, Qexc };   // Interaction types

	const G4int nCutMax;
	const G4double ThresholdParameter;
	const G4double QGSMThreshold;
	const G4double theNucleonRadius;

        G4ThreeVector theCurrentVelocity;                  // Uzhi 17 Apr. 2015 Is it needed?
        G4QGSMSplitableHadron* theProjectileSplitable;
private:
    G4ReactionProduct theProjectile;       

    G4Reggeons* Regge;     // Uzhi 18 Oct. 2016
    G4int InteractionMode;

    G4double alpha;
    G4double beta;

    G4double sigmaPt;   

    G4Nucleon* TheInvolvedNucleonsOfTarget[250];
    G4int NumberOfInvolvedNucleonsOfTarget;

    G4Nucleon* TheInvolvedNucleonsOfProjectile[250];
    G4int NumberOfInvolvedNucleonsOfProjectile;

    G4LorentzVector ProjectileResidual4Momentum;
    G4int           ProjectileResidualMassNumber;
    G4int           ProjectileResidualCharge;
    G4double        ProjectileResidualExcitationEnergy;

    G4LorentzVector TargetResidual4Momentum;
    G4int           TargetResidualMassNumber;
    G4int           TargetResidualCharge;
    G4double        TargetResidualExcitationEnergy;

private:
    // Parameters of nuclear destruction
    G4double CofNuclearDestruction;   // Cnd of nuclear destruction
    G4double R2ofNuclearDestruction;  // R2nd

    G4double ExcitationEnergyPerWoundedNucleon;

    G4double DofNuclearDestruction;       // D for momentum sampling
    G4double Pt2ofNuclearDestruction;     // Pt2
    G4double MaxPt2ofNuclearDestruction;  // Max Pt2

};

inline void G4QGSParticipants::StartPartonPairLoop()
{
}

inline G4PartonPair* G4QGSParticipants::GetNextPartonPair()
{
	if (thePartonPairs.empty()) return 0;
	G4PartonPair * result = thePartonPairs.back();
	thePartonPairs.pop_back();
	return result;
}


inline void G4QGSParticipants::SplitHadrons()
{
	unsigned int i;
	for(i = 0; i < theInteractions.size(); i++)
	{
		theInteractions[i]->SplitHadrons();
	}
}
//--------------------------------------
//Uzhi Copy from FTF Model.hh
/*
inline G4V3DNucleus* G4QGSParticipants::GetWoundedNucleus() const {
  return theNucleus;
}
*/

inline G4V3DNucleus* G4QGSParticipants::GetTargetNucleus() const {
  return theNucleus;
}

inline G4V3DNucleus* G4QGSParticipants::GetProjectileNucleus() const {
  return 0;
}

// Uzhi Start copy from FTFparameters
// Set parameters of nuclear destruction

inline void G4QGSParticipants::SetCofNuclearDestruction( const G4double aValue ) {
  CofNuclearDestruction = aValue;
}

inline void G4QGSParticipants::SetR2ofNuclearDestruction( const G4double aValue ) {
  R2ofNuclearDestruction = aValue;
}

inline void G4QGSParticipants::SetExcitationEnergyPerWoundedNucleon( const G4double aValue ) {
  ExcitationEnergyPerWoundedNucleon = aValue;
}

inline void G4QGSParticipants::SetDofNuclearDestruction( const G4double aValue ) {
  DofNuclearDestruction = aValue;
}

inline void G4QGSParticipants::SetPt2ofNuclearDestruction( const G4double aValue ) {
  Pt2ofNuclearDestruction = aValue;
}

inline void G4QGSParticipants::SetMaxPt2ofNuclearDestruction( const G4double aValue ) {
  MaxPt2ofNuclearDestruction = aValue;
}

// Get parameters of nuclear destruction
inline G4double G4QGSParticipants::GetCofNuclearDestruction() {
  return CofNuclearDestruction;
}

inline G4double G4QGSParticipants::GetR2ofNuclearDestruction() {
  return R2ofNuclearDestruction;
}

inline G4double G4QGSParticipants::GetExcitationEnergyPerWoundedNucleon() {
  return ExcitationEnergyPerWoundedNucleon;
}

inline G4double G4QGSParticipants::GetDofNuclearDestruction() {
  return DofNuclearDestruction;
}

inline G4double G4QGSParticipants::GetPt2ofNuclearDestruction() {
  return Pt2ofNuclearDestruction;
}

inline G4double G4QGSParticipants::GetMaxPt2ofNuclearDestruction() {
  return MaxPt2ofNuclearDestruction;
}
//Uzhi End copy from FTFparameters
#endif


