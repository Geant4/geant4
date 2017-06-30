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
// $Id: G4QGSDiffractiveExcitation.hh 102316 2017-01-20 16:12:52Z gcosmo $

#ifndef G4QGSDiffractiveExcitation_h
#define G4QGSDiffractiveExcitation_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QGSDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

// Modified:
//  25-05-07 : G.Folger
//       move from management/G4DiffractiveExcitation to to qgsm/G4QGSDiffractiveExcitation
//                  

#include "globals.hh"
class G4VSplitableHadron;
class G4ExcitedString;
#include "G4ThreeVector.hh"

class G4QGSDiffractiveExcitation 
{

public:

	G4QGSDiffractiveExcitation();
	virtual ~G4QGSDiffractiveExcitation();

	virtual G4bool ExciteParticipants (G4VSplitableHadron * aPartner, 
                                           G4VSplitableHadron * bPartner,        // Uzhi Oct. 2016
                                           G4bool ProjectileDiffraction=TRUE) const;  // Uzhi Oct. 2016   , G4bool ProjectileDiffraction
	virtual G4ExcitedString * String(G4VSplitableHadron * aHadron, G4bool isProjectile) const;

	//      void SetPtWidth(G4double aValue) { widthOfPtSquare = aValue*aValue; }
	//      void SetExtraMass(G4double aValue) { minExtraMass = aValue; }
	//      void SetMinimumMass(G4double aValue) { minmass = aValue; }


private:

	G4QGSDiffractiveExcitation(const G4QGSDiffractiveExcitation &right);

	G4double ChooseP(G4double Pmin, G4double Pmax) const;

	G4ThreeVector GaussianPt(G4double  AveragePt2, G4double maxPtSquare) const;

	const G4QGSDiffractiveExcitation & operator=(const G4QGSDiffractiveExcitation &right);
	int operator==(const G4QGSDiffractiveExcitation &right) const;
	int operator!=(const G4QGSDiffractiveExcitation &right) const;

private:
};

#endif
