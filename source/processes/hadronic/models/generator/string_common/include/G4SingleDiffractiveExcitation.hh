// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SingleDiffractiveExcitation.hh,v 1.2 1999/12/15 14:52:46 gunter Exp $

#ifndef G4SingleDiffractiveExcitation_h
#define G4SingleDiffractiveExcitation_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4SingleDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
class G4VSplitableHadron;
class G4ExcitedString;
#include "G4ThreeVector.hh"
#include "G4DiffractiveExcitation.hh"

class G4SingleDiffractiveExcitation : public G4DiffractiveExcitation
{

  public:

      G4SingleDiffractiveExcitation(G4double sigmaPt=0.6*GeV, G4double minExtraMass=250*MeV,G4double x0mass=250*MeV);
      ~G4SingleDiffractiveExcitation();

      G4bool ExciteParticipants (G4VSplitableHadron *aPartner, G4VSplitableHadron * bPartner) const;
      
//      void SetPtWidth(G4double aValue) { widthOfPtSquare = aValue*aValue; }
//      void SetExtraMass(G4double aValue) { minExtraMass = aValue; }
//      void SetMinimumMass(G4double aValue) { minmass = aValue; }


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
