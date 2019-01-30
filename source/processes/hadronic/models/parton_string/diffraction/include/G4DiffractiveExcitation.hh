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

#ifndef G4DiffractiveExcitation_h
#define G4DiffractiveExcitation_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//      Take a projectile and a target
//      excite the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"
#include "G4ThreeVector.hh"

class G4VSplitableHadron;
class G4ExcitedString;


class G4DiffractiveExcitation {
  public:
    G4DiffractiveExcitation();
    virtual ~G4DiffractiveExcitation();

    virtual G4bool ExciteParticipants( G4VSplitableHadron* aPartner, 
                                       G4VSplitableHadron* bPartner,
                                       G4FTFParameters* theParameters,
                                       G4ElasticHNScattering* theElastic ) const;

    virtual void CreateStrings( G4VSplitableHadron* aHadron, 
                                G4bool isProjectile,
                                G4ExcitedString*& FirstString, 
                                G4ExcitedString*& SecondString,
                                G4FTFParameters* theParameters ) const;

  private:
    G4DiffractiveExcitation( const G4DiffractiveExcitation& right );
    const G4DiffractiveExcitation& operator=( const G4DiffractiveExcitation& right );
    G4bool operator==( const G4DiffractiveExcitation& right ) const;
    G4bool operator!=( const G4DiffractiveExcitation& right ) const;

    G4double LambdaF(G4double sqrM, G4double sqrM1, G4double sqrM2) const;
      
    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;
    G4double ChooseP( G4double Pmin, G4double Pmax ) const;
    G4double GetQuarkFractionOfKink( G4double zmin, G4double zmax ) const;
    void UnpackMeson( G4int IdPDG, G4int& Q1, G4int& Q2 ) const;
    void UnpackBaryon( G4int IdPDG, G4int& Q1, G4int& Q2, G4int& Q3 ) const;
    G4int NewNucleonId( G4int Q1, G4int Q2, G4int Q3 ) const;

    // The "ExciteParticipants" method uses the following struct and 3 new utility methods:
    struct CommonVariables {
      G4int ProjectilePDGcode = 0, absProjectilePDGcode = 0, TargetPDGcode = 0, 
        absTargetPDGcode = 0;
      G4double M0projectile = 0.0, M0projectile2 = 0.0, M0target = 0.0, M0target2 = 0.0, 
        ProjMassT = 0.0, ProjMassT2 = 0.0, TargMassT = 0.0, TargMassT2 = 0.0, 
        MminProjectile = 0.0, MminTarget = 0.0, 
        ProjectileDiffStateMinMass = 0.0, ProjectileDiffStateMinMass2 = 0.0, 
        ProjectileNonDiffStateMinMass = 0.0, ProjectileNonDiffStateMinMass2 = 0.0,
        TargetDiffStateMinMass = 0.0, TargetDiffStateMinMass2 = 0.0, 
        TargetNonDiffStateMinMass = 0.0, TargetNonDiffStateMinMass2 = 0.0, 
        S = 0.0, SqrtS = 0.0, Pt2 = 0.0, PZcms = 0.0, PZcms2 = 0.0, 
        AveragePt2 = 0.0, maxPtSquare = 0.0,
        ProbExc = 0.0, Qminus = 0.0, Qplus = 0.0,
        PMinusNew = 0.0, PPlusNew = 0.0, TMinusNew = 0.0, TPlusNew = 0.0,
        PMinusMin = 0.0, PMinusMax = 0.0, TPlusMin = 0.0, TPlusMax = 0.0,
        ProbProjectileDiffraction = 0.0, ProbTargetDiffraction = 0.0, ProbOfDiffraction = 0.0;
      G4LorentzVector Pprojectile, Ptarget, Qmomentum;
      G4LorentzRotation toCms, toLab;
      G4SampleResonance BrW;
    };
    G4int ExciteParticipants_doChargeExchange( G4VSplitableHadron*    projectile,
                                               G4VSplitableHadron*    target,
                                               G4FTFParameters*       theParameters,
                                               G4ElasticHNScattering* theElastic,
                                               CommonVariables&       common ) const;
    G4bool ExciteParticipants_doDiffraction( G4VSplitableHadron* projectile,
                                             G4VSplitableHadron* target,
                                             G4FTFParameters*    theParameters,
                                             CommonVariables&    common ) const;
    G4bool ExciteParticipants_doNonDiffraction( G4VSplitableHadron* projectile,
                                                G4VSplitableHadron* target,
                                                G4FTFParameters*    theParameters,
                                                CommonVariables&    common ) const;
};

#endif

