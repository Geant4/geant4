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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//    Extends G4Track properties with information needed for the
//    channeling biasing operator.
//
// ------------------ G4ChannelingTrackData ------------------
//
// Author: M.Verderi (LLR), E. Bagli (INFN) April 2016
//
// --------------------------------------------------------------------

#ifndef G4ChannelingTrackData_hh
#define G4ChannelingTrackData_hh

class G4Channeling;
#include "G4VAuxiliaryTrackInformation.hh"
#include "G4ThreeVector.hh"

class G4ChannelingTrackData : public G4VAuxiliaryTrackInformation {
    friend class G4Channeling;
    
public:
    G4ChannelingTrackData();
    ~G4ChannelingTrackData();
    
    void Print() const;
    
private:
    const G4Channeling* fChannelingProcess;
    
public:
    void Reset(){
        fChannelingProcess = nullptr;
        fNuD = fElD = 1.;
        fPosCh = fMomCh = fDBL;
    }
    
public:
    G4double GetDensity() {return (fNuD + fElD) * 0.5;}
    
    void SetNuD(G4double aDouble) {fNuD = aDouble;};
    G4double GetNuD() {return fNuD;};
    
    void SetElD(G4double aDouble) {fElD = aDouble;};
    G4double GetElD() {return fElD;};

    void SetEFX(G4double aDouble) {fEFX = aDouble;};
    G4double GetEFX() {return fEFX;};

    void SetEFY(G4double aDouble) {fEFY = aDouble;};
    G4double GetEFY() {return fEFY;};
    
    G4ThreeVector GetMomCh() {return fMomCh;}
    void SetMomCh(G4ThreeVector a3vec) {fMomCh = a3vec;}
    
    G4ThreeVector GetPosCh() {return fPosCh;}
    void SetPosCh(G4ThreeVector a3vec) {fPosCh = a3vec;}
    
private:
    // ----------
    // DBL_MAX vector for reset
    // ----------
    G4ThreeVector fDBL;

    // ----------
    // Positiona and momentum in the crystal reference frame
    // ----------
    G4ThreeVector fMomCh;
    G4ThreeVector fPosCh;

    G4double fNuD;
    G4double fElD;

    G4double fEFX;
    G4double fEFY;
};

#endif
