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
// $Id: G4PreCompoundModel.hh,v 1.14 2002/01/22 13:47:28 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// by V. Lara

// Class Description
// Model implementation for pre-equilibrium decay models in geant4. 
// To be used in your physics list, in case you neeed this kind of physics.
// Can be used as a stand-allone model, but also in conjunction with an intra-nuclear
// transport, or any of the string-parton models.
// Class Description - End

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4PreCompoundEmission.hh"
#include "G4LorentzVector.hh"

#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"


//#define debug
//#define verbose
//#define pctest

class G4PreCompoundModel : public G4VPreCompoundModel
{
public:
    
    G4PreCompoundModel(G4ExcitationHandler * const value) : 
	G4VPreCompoundModel(value) {};

    ~G4PreCompoundModel() {};

private:
    G4PreCompoundModel() {};

    G4PreCompoundModel(const G4PreCompoundModel &right) {};

    const G4PreCompoundModel& operator=(const G4PreCompoundModel &right);
    G4bool operator==(const G4PreCompoundModel &right) const;
    G4bool operator!=(const G4PreCompoundModel &right) const;

public:
    G4VParticleChange * ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus);

    G4ReactionProductVector* DeExcite(const G4Fragment& aFragment) const;

#ifdef pctest
    static G4Fragment GetInitialFragmentForTest()
	{ return G4PreCompoundModel::theInitialFragmentForTest; }
    static G4std::vector<G4String*> * GetCreatorModels()
	{ return &G4PreCompoundModel::theCreatorModels; }
#endif

private:  

    void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				    G4ReactionProductVector * theResult) const;

#ifdef debug				  
    void CheckConservation(const G4Fragment & theInitialState,
			   const G4Fragment & aFragment,
			   G4ReactionProductVector * Result) const;
#endif

    G4ParticleChange theResult;

#ifdef pctest
    static G4Fragment theInitialFragmentForTest;
    static G4std::vector<G4String*> theCreatorModels;
#endif

};
#endif











