//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DiffractiveSplitableHadron.hh,v 1.5 2001/08/01 17:05:44 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
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

class G4DiffractiveSplitableHadron : public G4VSplitableHadron
{

public:

	G4DiffractiveSplitableHadron(const G4ReactionProduct & aPrimary);
	G4DiffractiveSplitableHadron(const G4Nucleon & aNucleon);
	G4DiffractiveSplitableHadron(const G4VKineticNucleon * aNucleon);
	~G4DiffractiveSplitableHadron();


	int operator==(const G4DiffractiveSplitableHadron &right) const;
	int operator!=(const G4DiffractiveSplitableHadron &right) const;


	void SplitUp();
	G4Parton * GetNextParton() ;
	G4Parton * GetNextAntiParton();
	
private:
	G4DiffractiveSplitableHadron();
	G4DiffractiveSplitableHadron(const G4DiffractiveSplitableHadron &right);
	const G4DiffractiveSplitableHadron & operator=(const G4DiffractiveSplitableHadron &right);

//implementation
	G4int Diquark(G4int aquark,G4int bquark,G4int Spin) const; // to splitable hadron
	void ChooseStringEnds(G4int PDGcode,G4int * aEnd, G4int * bEnd) const; // to splitable hadron

private:
	G4Parton *Parton[2];
	G4int    PartonIndex; 

};

#endif	
