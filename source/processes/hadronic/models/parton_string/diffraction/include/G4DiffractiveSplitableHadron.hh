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
// $Id: G4DiffractiveSplitableHadron.hh,v 1.6 2010-09-20 15:50:46 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
