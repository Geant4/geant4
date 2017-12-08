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
// $Id: G4SingleDiffractiveExcitation.hh 106980 2017-10-31 09:02:49Z gcosmo $

#ifndef G4SingleDiffractiveExcitation_h
#define G4SingleDiffractiveExcitation_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4SingleDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4QGSDiffractiveExcitation.hh"

class G4VSplitableHadron;
class G4ExcitedString;

class G4SingleDiffractiveExcitation : public G4QGSDiffractiveExcitation
{
  public:

    G4SingleDiffractiveExcitation(G4double sigmaPt=0.6*CLHEP::GeV,
                                  G4double minExtraMass=250*CLHEP::MeV,
                                  G4double x0mass=250*CLHEP::MeV);
    ~G4SingleDiffractiveExcitation();

    G4bool ExciteParticipants(G4VSplitableHadron *aPartner, G4VSplitableHadron * bPartner) const;

    //void SetPtWidth(G4double aValue) { widthOfPtSquare = aValue*aValue; }
    //void SetExtraMass(G4double aValue) { minExtraMass = aValue; }
    //void SetMinimumMass(G4double aValue) { minmass = aValue; }

  private:

    G4SingleDiffractiveExcitation(const G4SingleDiffractiveExcitation &right);

    G4double ChooseX(G4double Xmin, G4double Xmax) const;
    G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;

    const G4SingleDiffractiveExcitation & operator=(const G4SingleDiffractiveExcitation &right);
    int operator==(const G4SingleDiffractiveExcitation &right) const;
    int operator!=(const G4SingleDiffractiveExcitation &right) const;

  private:
    // Model Parameters:
    const G4double widthOfPtSquare;	// width^2 of pt for string excitation
    const G4double minExtraMass;	// minimum excitation mass 
    const G4double minmass;	// mean pion transverse mass; used for Xmin 
};

#endif

