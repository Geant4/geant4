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
// $Id: G4DiffractiveSplitableHadron.hh 100828 2016-11-02 15:25:59Z gcosmo $
// GEANT4 tag $Name:  $
//

#ifndef G4DiffractiveSplitableHadron_h
#define G4DiffractiveSplitableHadron_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4DiffractiveSplitableHadron----------------
//             by Gunter Folger, August 1998.
//       class splitting an interacting particle. Used by FTF String Model.
// ------------------------------------------------------------

#include "G4VSplitableHadron.hh"
#include "G4Nucleon.hh"
#include "G4Parton.hh"


class G4DiffractiveSplitableHadron : public G4VSplitableHadron {
  public:
    G4DiffractiveSplitableHadron();
    G4DiffractiveSplitableHadron( const G4ReactionProduct& aPrimary );
    G4DiffractiveSplitableHadron( const G4Nucleon& aNucleon );
    G4DiffractiveSplitableHadron( const G4VKineticNucleon* aNucleon );
    ~G4DiffractiveSplitableHadron();

    void SplitUp();
    G4Parton* GetNextParton() ;
    G4Parton* GetNextAntiParton();

    void SetFirstParton( G4int PDGcode );
    void SetSecondParton( G4int PDGcode );

  private:
    G4DiffractiveSplitableHadron( const G4DiffractiveSplitableHadron& );
    G4DiffractiveSplitableHadron& operator=( const G4DiffractiveSplitableHadron& );
    int operator==( const G4DiffractiveSplitableHadron& right ) const;
    int operator!=( const G4DiffractiveSplitableHadron& right ) const;

    G4int Diquark( G4int aquark, G4int bquark, G4int Spin ) const;
    void ChooseStringEnds( G4int PDGcode, G4int* aEnd, G4int* bEnd ) const;

    G4Parton* Parton[2];
    G4int PartonIndex; 
};

#endif

