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
#ifndef G4QGSMSplitableHadron_h
#define G4QGSMSplitableHadron_h 1

#include "G4VSplitableHadron.hh"
#include "G4PartonVector.hh"
#include "G4MesonSplitter.hh"
#include "G4BaryonSplitter.hh"
#include "Randomize.hh"
#include <deque>

// based on prototype by Maxim Komogorov
// Splitting into methods, and centralizing of model parameters HPW Feb 1999
// continued clean-up of interfaces and algorithms HPW 1999.
// Redesign of data structures and algorithms HPW Feb 1999

class G4QGSMSplitableHadron : public G4VSplitableHadron
{

public:
	G4QGSMSplitableHadron();
	G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary);
	G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary, G4bool Direction);
	G4QGSMSplitableHadron(const G4Nucleon & aNucleon);
	G4QGSMSplitableHadron(const G4Nucleon & aNucleon, G4bool Direction);

	virtual ~G4QGSMSplitableHadron();

private:
	const G4QGSMSplitableHadron & operator=(const G4QGSMSplitableHadron &right);

public:
	virtual void SplitUp();
	virtual void SetFirstParton(G4int PDGcode);  // Uzhi 24.11.10
	virtual void SetSecondParton(G4int PDGcode);  // Uzhi 24.11.10
	virtual G4Parton * GetNextParton();
	virtual G4Parton * GetNextAntiParton();

private:
	void InitParameters();
	void DiffractiveSplitUp();
	void SoftSplitUp();

	G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare);
	void GetValenceQuarkFlavors(const G4ParticleDefinition * aPart,
			G4Parton *& Parton1, G4Parton *& Parton2);
	G4Parton * BuildSeaQuark(G4bool isAntiQuark, G4int aPDGCode, G4int nSeaPair);
	G4double SampleX(G4double anXmin, G4int nSea, G4int theTotalSea, G4double aBeta);

private:
	// aggregated data
	G4bool Direction; // FALSE is target. - candidate for more detailed design. @@@@ HPW

	std::deque<G4Parton *> Color;
	std::deque<G4Parton *> AntiColor;
//std::deque<G4Parton *>::iterator iP;      // Uzhi
//std::deque<G4Parton *>::iterator iAP;     // Uzhi
unsigned int iP;                            // Uzhi 5.06.2015
unsigned int iAP;                           // Uzhi 5.06.2015
private:

	// associated classes
	G4MesonSplitter theMesonSplitter;
	G4BaryonSplitter theBaryonSplitter;

private:
	// model parameters
	G4double alpha;
	G4double beta;
	G4double theMinPz;
	G4double StrangeSuppress;
	G4double sigmaPt;
	G4double widthOfPtSquare;
	G4double minTransverseMass;
};

inline G4Parton* G4QGSMSplitableHadron::GetNextParton()
{
	if(Color.size()==0) return 0;

        G4Parton * result = Color.operator[](iP);
        iP++; if(iP == Color.size()) iP=0;
	return result;
}

inline G4Parton* G4QGSMSplitableHadron::GetNextAntiParton()
{
	if(AntiColor.size() == 0) return 0;

        G4Parton * result = AntiColor.operator[](iAP);
        iAP++; if(iAP == AntiColor.size()) iAP=0;
	return result;
}

inline void G4QGSMSplitableHadron::SetFirstParton(G4int PDGcode)
{PDGcode++;}
inline void G4QGSMSplitableHadron::SetSecondParton(G4int PDGcode)
{PDGcode++;}
#endif


