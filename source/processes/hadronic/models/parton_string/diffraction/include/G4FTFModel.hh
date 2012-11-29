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
// $Id$
// GEANT4 tag $Name:  $
//
// Class Description
// Final state production code for hadron inelastic scattering above 20 GeV
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

class G4VSplitableHadron;
class G4ExcitedString;

#include "G4FTFParameters.hh"
#include "G4FTFParticipants.hh"

#include "G4ExcitedStringVector.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4ElasticHNScattering.hh"
#include "G4FTFAnnihilation.hh"

class G4FTFModel : public G4VPartonStringModel
{

  public:
      G4FTFModel(const G4String& modelName = "FTF");
      //G4FTFModel(G4double , G4double , G4double );
      //G4FTFModel(G4DiffractiveExcitation * anExcitation);
      ~G4FTFModel();
  private:
      G4FTFModel(const G4FTFModel &right);
      const G4FTFModel & operator=(const G4FTFModel &right);
      int operator==(const G4FTFModel &right) const;
      int operator!=(const G4FTFModel &right) const;
  public:
      void Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile);
      G4ExcitedStringVector * GetStrings();
      G4V3DNucleus * GetWoundedNucleus() const;
      virtual void ModelDescription(std::ostream&) const;

  protected:
  
  private:
       void StoreInvolvedNucleon();              
       void ReggeonCascade();
       G4bool PutOnMassShell();
       G4bool ExciteParticipants();
       G4ExcitedStringVector * BuildStrings();
       void GetResidualNucleus();                  
       void AjustTargetNucleonForAnnihilation(G4VSplitableHadron *SelectedAntiBaryon,
                                              G4VSplitableHadron *SelectedTargetNucleon); 
       G4ThreeVector GaussianPt(G4double  AveragePt2, G4double maxPtSquare) const;
  
  private:     

       G4ReactionProduct theProjectile;       
       G4FTFParticipants theParticipants;
       
       G4Nucleon * TheInvolvedNucleon[250];
       G4int  NumberOfInvolvedNucleon;
       G4int  NumberOfInvolvedTargetNucleon;

       G4Nucleon * TheInvolvedNucleonOfProjectile[250];
       G4int  NumberOfInvolvedNucleonOfProjectile;

       G4FTFParameters  *theParameters;
       G4DiffractiveExcitation * theExcitation;
       G4ElasticHNScattering   * theElastic;
       G4FTFAnnihilation       * theAnnihilation;              // Uzhi 17.11.10

       std::vector<G4VSplitableHadron *> theAdditionalString;  // Uzhi 17.11.10

       G4LorentzVector Residual4Momentum;
       G4double ResidualExcitationEnergy;

};

// ------------------------------------------------------------
inline 
G4V3DNucleus * G4FTFModel::GetWoundedNucleus() const
{
	return theParticipants.GetWoundedNucleus();
}

#endif
