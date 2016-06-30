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
#include "G4PomeronCrossSection.hh"
#include "G4QGSDiffractiveExcitation.hh"
#include "G4SingleDiffractiveExcitation.hh"
#include "G4PartonPair.hh" 
#include "G4QGSMSplitableHadron.hh" 
#include "G4V3DNucleus.hh"

#include "G4VSplitableHadron.hh"  // Uzhi

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
                theCurrentVelocity = -aBoost;                      // Uzhi 17 Apr. 2015
		if(theNucleus) theNucleus->DoLorentzBoost(aBoost);
		theBoost = aBoost;
	}

	G4PartonPair* GetNextPartonPair();
	void BuildInteractions(const G4ReactionProduct  &thePrimary);
	void StartPartonPairLoop();

//Uzhi Start copy from FTFmodel
private:
//Uzhi    G4V3DNucleus* GetWoundedNucleus() const;
    G4V3DNucleus* GetTargetNucleus() const;
    G4V3DNucleus* GetProjectileNucleus() const;

    void PrepareInitialState( const G4ReactionProduct& thePrimary );
    void GetList( const G4ReactionProduct& thePrimary );

    void StoreInvolvedNucleon();              
    void ReggeonCascade();
    G4bool PutOnMassShell();
//Uzhi    G4bool ExciteParticipants();
//Uzhi    G4ExcitedStringVector* BuildStrings();
    void GetResiduals();                  

//Uzhi    G4bool AdjustNucleons( G4VSplitableHadron* SelectedAntiBaryon,
//Uzhi                           G4Nucleon*          ProjectileNucleon,
//Uzhi                           G4VSplitableHadron* SelectedTargetNucleon,
//Uzhi                           G4Nucleon*          TargetNucleon,
//Uzhi                           G4bool              Annihilation ); 
    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;

    G4bool ComputeNucleusProperties( G4V3DNucleus* nucleus, G4LorentzVector& nucleusMomentum, 
                                     G4LorentzVector& residualMomentum, G4double& sumMasses,   
                                     G4double& residualExcitationEnergy, G4double& residualMass,
                                     G4int& residualMassNumber, G4int& residualCharge );
    // Utility method used by PutOnMassShell.

    G4bool GenerateDeltaIsobar( const G4double sqrtS, const G4int numberOfInvolvedNucleons,
                                G4Nucleon* involvedNucleons[], G4double& sumMasses );
    // Utility method used by PutOnMassShell.

    G4bool SamplingNucleonKinematics( G4double averagePt2, const G4double maxPt2,
                                      G4double dCor, G4V3DNucleus* nucleus, 
                                      const G4LorentzVector& pResidual, 
                                      const G4double residualMass, const G4int residualMassNumber,
                                      const G4int numberOfInvolvedNucleons,
                                      G4Nucleon* involvedNucleons[], G4double& mass2 );
    // Utility method used by PutOnMassShell.

    G4bool CheckKinematics( const G4double sValue, const G4double sqrtS, 
                            const G4double projectileMass2, const G4double targetMass2,
                            const G4double nucleusY, const G4bool isProjectileNucleus,
                            const G4int numberOfInvolvedNucleons, G4Nucleon* involvedNucleons[],
                            G4double& targetWminus, G4double& projectileWplus, G4bool& success );
    // Utility method used by PutOnMassShell.

    G4bool FinalizeKinematics( const G4double w, const G4bool isProjectileNucleus, 
                               const G4LorentzRotation& boostFromCmsToLab,
                               const G4double residualMass, const G4int residualMassNumber,
                               const G4int numberOfInvolvedNucleons, 
                               G4Nucleon* involvedNucleons[],
	                       G4LorentzVector& residual4Momentum );
    // Utility method used by PutOnMassShell.

//Uzhi End copy from FTFmodel

void CreateStrings();

//Uzhi Start copy from FTFparameters
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
//Uzhi End copy from FTFparameters

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

	G4SingleDiffractiveExcitation theSingleDiffExcitation;
	G4QGSDiffractiveExcitation theDiffExcitaton;
	G4int ModelMode;
	G4bool IsSingleDiffractive();

	G4ThreeVector theBoost;
	G4double SampleX(G4double anXmin, G4int nSea, G4int theTotalSea, G4double aBeta);
protected:
	// model parameters HPW
	enum  { SOFT, DIFFRACTIVE };
	const G4int nCutMax;
	const G4double ThresholdParameter;
	const G4double QGSMThreshold;
	const G4double theNucleonRadius;

   // cash theCurrentVelocity for lorentztrafo HPW 
        G4ThreeVector theCurrentVelocity;                  // Uzhi 17 Apr. 2015
    G4QGSMSplitableHadron* theProjectileSplitable;         // Uzhi 21.05.2015
private:
//Uzhi Start copy from FTFmodel
    G4ReactionProduct theProjectile;       
//    G4QGSMSplitableHadron* theProjectileSplitable;

    G4double alpha;
    G4double beta;

    G4double sigmaPt;
//Uzhi    G4FTFParticipants theParticipants;       

    G4Nucleon* TheInvolvedNucleonsOfTarget[250];
    G4int NumberOfInvolvedNucleonsOfTarget;

    G4Nucleon* TheInvolvedNucleonsOfProjectile[250];
    G4int NumberOfInvolvedNucleonsOfProjectile;

//Uzhi    G4FTFParameters* theParameters;
//Uzhi    G4DiffractiveExcitation* theExcitation;
//Uzhi    G4ElasticHNScattering* theElastic;
//Uzhi    G4FTFAnnihilation* theAnnihilation;  

//Uzhi    std::vector< G4VSplitableHadron* > theAdditionalString; 

//Uzhi    G4double LowEnergyLimit;
//Uzhi    G4bool HighEnergyInter;

    G4LorentzVector ProjectileResidual4Momentum;
    G4int           ProjectileResidualMassNumber;
    G4int           ProjectileResidualCharge;
    G4double        ProjectileResidualExcitationEnergy;

    G4LorentzVector TargetResidual4Momentum;
    G4int           TargetResidualMassNumber;
    G4int           TargetResidualCharge;
    G4double        TargetResidualExcitationEnergy;
//Uzhi End copy from FTFmodel

//Uzhi Start copy from FTFparameters
private:
    // Parameters of nuclear destruction
    G4double CofNuclearDestruction;   // Cnd of nuclear destruction
    G4double R2ofNuclearDestruction;  // R2nd

    G4double ExcitationEnergyPerWoundedNucleon;

    G4double DofNuclearDestruction;       // D for momentum sampling
    G4double Pt2ofNuclearDestruction;     // Pt2
    G4double MaxPt2ofNuclearDestruction;  // Max Pt2
//Uzhi End copy from FTFparameters

};


inline G4bool G4QGSParticipants::IsSingleDiffractive()
{
	G4bool result=false;
	if(G4UniformRand()<1.) result = true;
	return result;
}

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
//G4cout<<"----------------------------------- SplitHadrons -------------"<<G4endl;
//G4cout<<"theInteractions.size() "<<theInteractions.size()<<G4endl;

	unsigned int i;
	for(i = 0; i < theInteractions.size(); i++)
	{
//G4cout<<i<<" "<<theInteractions[i]->GetProjectile()<<" "<<theInteractions[i]->GetProjectileNucleon()<<" "
//              <<theInteractions[i]->GetTarget()<<" "<<theInteractions[i]->GetTargetNucleon()<<G4endl;

//G4cout<<i<<" "<<theInteractions[i]->GetProjectile()->GetDefinition()->GetParticleName()<<G4endl;
//G4cout<<i<<" "<<theInteractions[i]->GetTarget()->GetDefinition()->GetParticleName()<<G4endl;
	}
	for(i = 0; i < theInteractions.size(); i++)
	{
//G4cout<<i<<" "<<theInteractions[i]->GetNumberOfSoftCollisions()<<" "
//              <<theInteractions[i]->GetNumberOfHardCollisions()<<" "
//              <<theInteractions[i]->GetNumberOfDiffractiveCollisions()<<G4endl;

	}
//G4cout<<"******************************** SplitHadrons ************"<<G4endl;
	for(i = 0; i < theInteractions.size(); i++)
	{
//G4cout<<"Interaction # "<<i<<" QGSPartic"<<G4endl;
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

//Uzhi Start copy from FTFparameters


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


