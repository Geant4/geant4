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
// $Id: G4DiffractiveExcitation.hh,v 1.3 2007-04-24 10:37:10 gunter Exp $

#ifndef G4DiffractiveExcitation_h
#define G4DiffractiveExcitation_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
class G4VSplitableHadron;
class G4ExcitedString;
#include "G4ThreeVector.hh"

class G4DiffractiveExcitation 
{

  public:

      G4DiffractiveExcitation();                           // Uzhi
      virtual ~G4DiffractiveExcitation();

      virtual G4bool ExciteParticipants (G4VSplitableHadron *aPartner, G4VSplitableHadron * bPartner) const;
      virtual G4ExcitedString * String(G4VSplitableHadron * aHadron, G4bool isProjectile) const;
      
//      void SetPtWidth(G4double aValue) { widthOfPtSquare = aValue*aValue; }
//      void SetExtraMass(G4double aValue) { minExtraMass = aValue; }
//      void SetMinimumMass(G4double aValue) { minmass = aValue; }


  private:

      G4DiffractiveExcitation(const G4DiffractiveExcitation &right);
      
//    G4double ChooseX(G4double Xmin, G4double Xmax) const;                       // Uzhi
      G4double ChooseP(G4double Pmin, G4double Pmax) const;                       // Uzhi

//    G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;
      G4ThreeVector GaussianPt(G4double  AveragePt2, G4double maxPtSquare) const;  // Uzhi
      
      const G4DiffractiveExcitation & operator=(const G4DiffractiveExcitation &right);
      int operator==(const G4DiffractiveExcitation &right) const;
      int operator!=(const G4DiffractiveExcitation &right) const;

  private:
// Model Parameters:
/*                                                                                // Uzhi
	const G4double widthOfPtSquare;	// width^2 of pt for string excitation    
	const G4double minExtraMass;	// minimum excitation mass                
	const G4double minmass;	// mean pion transverse mass; used for Xmin       
*/                                                                                // Uzhi
};

#endif
