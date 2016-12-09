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
// $Id: G4DiffractiveExcitation.hh 100828 2016-11-02 15:25:59Z gcosmo $

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
#include "G4ElasticHNScattering.hh"  // Uzhi 3.09.09
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
    int operator==( const G4DiffractiveExcitation& right ) const;
    int operator!=( const G4DiffractiveExcitation& right ) const;

    G4double LambdaF(G4double sqrM, G4double sqrM1, G4double sqrM2) const; // May 2014
      
    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;
    G4double ChooseP( G4double Pmin, G4double Pmax ) const;
    G4double GetQuarkFractionOfKink( G4double zmin, G4double zmax ) const;
    void UnpackMeson( G4int IdPDG, G4int& Q1, G4int& Q2 ) const;
    void UnpackBaryon( G4int IdPDG, G4int& Q1, G4int& Q2, G4int& Q3 ) const;
    G4int NewNucleonId( G4int Q1, G4int Q2, G4int Q3 ) const;
};

#endif

