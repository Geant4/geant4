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
// $Id: G4QuarkExchange.hh 66241 2012-12-13 18:34:42Z gunter $

#ifndef G4QuarkExchange_h
#define G4QuarkExchange_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QuarkExchange --------------
//             by V. Uzhinsky, October 2016.
//       QuarkExchange is used by strings models.
//	    Take a projectile and a target.
//Simulate Q exchange with excitation of projectile or target.
// ------------------------------------------------------------

//#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
//#include "G4QGSDiffractiveExcitation.hh"

class G4VSplitableHadron;
class G4ExcitedString;

class G4QuarkExchange
{

public:

	G4QuarkExchange();
	~G4QuarkExchange();

	G4bool ExciteParticipants (G4VSplitableHadron *aPartner, G4VSplitableHadron * bPartner) const;

private:
	G4QuarkExchange(const G4QuarkExchange &right);

	G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;

	const G4QuarkExchange & operator=(const G4QuarkExchange &right);
	int operator==(const G4QuarkExchange &right) const;
	int operator!=(const G4QuarkExchange &right) const;

};

#endif
