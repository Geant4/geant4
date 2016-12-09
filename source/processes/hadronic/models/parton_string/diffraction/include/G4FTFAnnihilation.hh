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
// $Id: G4FTFAnnihilation.hh 100828 2016-11-02 15:25:59Z gcosmo $

#ifndef G4FTFAnnihilation_h
#define G4FTFAnnihilation_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4FTFAnnihilation       --------------
//             by Vladimir Uzhinsky, November 2010.
//          Annihilation used by Fritiof (FTF) model
//              Takes a projectile and a target
//              Produces strings (excited hadrons)
// ------------------------------------------------------------

#include "globals.hh"
#include "G4FTFParameters.hh"
#include "G4ThreeVector.hh"

class G4VSplitableHadron;
class G4ExcitedString;


class G4FTFAnnihilation {
  public:
    G4FTFAnnihilation();
    virtual ~G4FTFAnnihilation();

    virtual G4bool Annihilate( G4VSplitableHadron* aPartner, 
                               G4VSplitableHadron* bPartner,
                               G4VSplitableHadron*& AdditionalString, 
                               G4FTFParameters* theParameters ) const;

  private:
    G4FTFAnnihilation( const G4FTFAnnihilation& right );      
    const G4FTFAnnihilation& operator=( const G4FTFAnnihilation& right );
    int operator==( const G4FTFAnnihilation& right ) const;
    int operator!=( const G4FTFAnnihilation& right ) const;

    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;
    G4double ChooseX( G4double Alpha, G4double Beta ) const;
    void UnpackBaryon( G4int IdPDG, G4int& Q1, G4int& Q2, G4int& Q3 ) const;
};

#endif

