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
// $Id: G4FTFModel.hh 100828 2016-11-02 15:25:59Z gcosmo $
// GEANT4 tag $Name:  $
//
// Class Description
// Final state production code for hadron inelastic scattering above 3 GeV
// based on the modeling ansatz used in FRITIOF.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with an object
// of G4TheoFSGenerator. 
// Class Description - End

#ifndef G4FTFModel_h
#define G4FTFModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------

#include "G4VPartonStringModel.hh"
#include "G4FTFParameters.hh"
#include "G4FTFParticipants.hh"
#include "G4ExcitedStringVector.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4ElasticHNScattering.hh"
#include "G4FTFAnnihilation.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

class G4VSplitableHadron;
class G4ExcitedString;


class G4FTFModel : public G4VPartonStringModel {
  public:
    G4FTFModel( const G4String& modelName = "FTF" );
    ~G4FTFModel();

    void Init( const G4Nucleus& aNucleus, const G4DynamicParticle& aProjectile );
    G4ExcitedStringVector* GetStrings();
    G4V3DNucleus* GetWoundedNucleus() const;
    G4V3DNucleus* GetTargetNucleus() const;
    G4V3DNucleus* GetProjectileNucleus() const;

    virtual void ModelDescription( std::ostream& ) const;

  private:
    G4FTFModel( const G4FTFModel& right );
    const G4FTFModel& operator=( const G4FTFModel& right );
    int operator==( const G4FTFModel& right ) const;
    int operator!=( const G4FTFModel& right ) const;

    void StoreInvolvedNucleon();              
    void ReggeonCascade();
    G4bool PutOnMassShell();
    G4bool ExciteParticipants();
    G4ExcitedStringVector* BuildStrings();
    void GetResiduals();                  
    G4bool AdjustNucleons( G4VSplitableHadron* SelectedAntiBaryon,
                           G4Nucleon*          ProjectileNucleon,
                           G4VSplitableHadron* SelectedTargetNucleon,
                           G4Nucleon*          TargetNucleon,
                           G4bool              Annihilation ); 
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

    G4ReactionProduct theProjectile;       
    G4FTFParticipants theParticipants;
       
    G4Nucleon* TheInvolvedNucleonsOfTarget[250];
    G4int NumberOfInvolvedNucleonsOfTarget;

    G4Nucleon* TheInvolvedNucleonsOfProjectile[250];
    G4int NumberOfInvolvedNucleonsOfProjectile;

    G4FTFParameters* theParameters;
    G4DiffractiveExcitation* theExcitation;
    G4ElasticHNScattering* theElastic;
    G4FTFAnnihilation* theAnnihilation;  

    std::vector< G4VSplitableHadron* > theAdditionalString; 

    G4double LowEnergyLimit;
    G4bool HighEnergyInter;

    G4LorentzVector ProjectileResidual4Momentum;
    G4int           ProjectileResidualMassNumber;
    G4int           ProjectileResidualCharge;
    G4double        ProjectileResidualExcitationEnergy;

    G4LorentzVector TargetResidual4Momentum;
    G4int           TargetResidualMassNumber;
    G4int           TargetResidualCharge;
    G4double        TargetResidualExcitationEnergy;
};


inline G4V3DNucleus* G4FTFModel::GetWoundedNucleus() const {
  return theParticipants.GetWoundedNucleus();
}


inline G4V3DNucleus* G4FTFModel::GetTargetNucleus() const {
  return theParticipants.GetWoundedNucleus();
}


inline G4V3DNucleus* G4FTFModel::GetProjectileNucleus() const {
  return theParticipants.GetProjectileNucleus();
}

#endif

