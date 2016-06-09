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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4InelasticSplitableHadron.hh,v 1.1 2003/10/07 11:26:25 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//

#ifndef G4InelasticSplitableHadron_h
#define G4InelasticSplitableHadron_h 1

#include "G4VSplitableHadron.hh"

class G4InelasticSplitableHadron : public G4VSplitableHadron
{

public:

	G4InelasticSplitableHadron(const G4ReactionProduct & aPrimary);
	G4InelasticSplitableHadron(const G4Nucleon & aNucleon);

	~G4InelasticSplitableHadron();

	int operator==(const G4InelasticSplitableHadron &right) const;
	int operator!=(const G4InelasticSplitableHadron &right) const;


	void SplitUp();
	G4Parton * GetNextParton() ;
	
private:
	G4InelasticSplitableHadron();
	G4InelasticSplitableHadron(const G4InelasticSplitableHadron &right);
	const G4InelasticSplitableHadron & operator=(const G4InelasticSplitableHadron &right);

//implementation
        void GetValenceQuarkFlavors(G4int PDGcode, G4int& aEnd, G4int& bEnd);
	G4int Diquark(G4int aquark,G4int bquark,G4int Spin) const; // to splitable hadron

private:
	G4Parton *Parton[2];
	G4int    PartonIndex; 

    G4double udspin1;         // 1./6.
    G4double uuspin1;         // 1./3.
    G4double udspin0;         // 1./2.

};

#endif	
